"""Declared external-cohort evaluation; evidence gaps remain first-class results."""
from __future__ import annotations
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from microglia_repurposing.contracts import (CLAIM_BOUNDARY, DISEASE_SCHEMA,
    file_sha256, require, sha256_identity, validate_gene_row)
from microglia_repurposing.identity import pin_symbol


def verify_document(document):
    require(document.get('content_hash') == sha256_identity(
        {k:v for k,v in document.items() if k != 'content_hash'}), 'document_hash_mismatch')
    return document


def seal(document):
    document['content_hash'] = sha256_identity({k:v for k,v in document.items() if k != 'content_hash'})
    return document


def validate_signature(bundle):
    verify_document(bundle)
    require(bundle.get('schema') == DISEASE_SCHEMA and bundle.get('research_planning_only') is True,
            'external_signature_contract_invalid')
    require(bundle.get('species') == 'human', 'external_signature_species_invalid')
    require(bool(bundle.get('source', {}).get('sha256')), 'external_signature_source_missing')
    seen = set()
    for gene in bundle['genes']:
        validate_gene_row(gene)
        identity = pin_symbol(gene['gene_symbol'])
        require(identity is not None and identity['hgnc_id'] == gene['hgnc_id'], 'external_gene_identity_mismatch')
        require(gene['hgnc_id'] not in seen, 'external_gene_duplicate')
        seen.add(gene['hgnc_id'])
        weight = gene['weight']
        require(isinstance(weight, (float, int)) and math.isfinite(weight) and weight != 0,
                'external_weight_invalid')
        require(gene['sign'] == (1 if weight > 0 else -1), 'external_direction_invalid')
    require(sum(g['sign'] == 1 for g in bundle['genes']) >= 3 and
            sum(g['sign'] == -1 for g in bundle['genes']) >= 3, 'external_signature_too_small')


def bh_adjust(values):
    """BH across the declared finite p-value family, preserving input order."""
    p = np.asarray(values, dtype=float)
    require(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(), 'invalid_p_value')
    if not len(p):
        return p
    order = np.argsort(p, kind='stable')
    ranked = np.minimum.accumulate((p[order] * len(p) / np.arange(1,len(p)+1))[::-1])[::-1]
    adjusted = np.empty(len(p)); adjusted[order] = np.minimum(ranked, 1)
    return adjusted


def seaad_signature(path: Path, state: str, protocol: dict):
    verify_document(protocol)
    require(state in protocol['states'], 'undeclared_state')
    receipt = json.loads(path.with_suffix('.provenance.json').read_text())
    require(file_sha256(path) == receipt['sha256'], 'external_source_hash_mismatch')
    require(receipt['protocol_hash'] == protocol['content_hash'], 'source_protocol_mismatch')
    table = pd.read_csv(path)
    effect = 'logFC_' + protocol['contrast']; pcolumn = 'p_' + protocol['contrast']
    require({'Unnamed: 0', 'gene_id', effect, pcolumn}.issubset(table.columns), 'external_columns_missing')
    require(not table['gene_id'].duplicated().any(), 'external_gene_id_duplicate')
    p = pd.to_numeric(table[pcolumn], errors='raise')
    effects = pd.to_numeric(table[effect], errors='raise')
    require(not np.isinf(p).any() and not np.isinf(effects).any(), 'infinite_external_value')
    finite = p.notna()
    table['q'] = np.nan
    table.loc[finite, 'q'] = bh_adjust(p[finite])
    table['effect'] = effects
    candidates = table[finite & effects.notna() & (effects != 0) & (table['q'] <= protocol['selection']['q_max'])].copy()
    candidates['magnitude'] = candidates['effect'].abs()
    candidates = candidates.sort_values(['q','magnitude','gene_id'], ascending=[True,False,True])
    genes = []; unmapped = []; seen = set(); counts = {1:0,-1:0}
    # Identity resolution precedes the frozen per-direction cap. No collapse of aliases.
    for _, row in candidates.iterrows():
        # Publisher gene_id is a numeric row index; symbols live in the unnamed CSV index.
        hit = pin_symbol(str(row['Unnamed: 0']))
        if hit is None:
            unmapped.append(str(row['Unnamed: 0'])); continue
        require(hit['hgnc_id'] not in seen, 'external_gene_duplicate')
        seen.add(hit['hgnc_id'])
        weight = float(row['effect']); sign = 1 if weight > 0 else -1
        if counts[sign] >= protocol['selection']['max_genes_per_direction']:
            continue
        counts[sign] += 1
        genes.append({**hit,'weight':weight,'sign':sign,'q':float(row['q']),'spine_entity_id':None})
    bundle = seal(dict(schema=DISEASE_SCHEMA, research_planning_only=True, claim_boundary=CLAIM_BOUNDARY,
        species='human', cell_context='human_brain_microglia', cell_state=state + '_pathology_pseudoprogression',
        source={'library':'SEA-AD MTG NEBULA published aggregate', 'url':receipt['url'],
                'sha256':receipt['sha256'],'doi':'10.1038/s41593-024-01774-5',
                'model':'published_NEBULA',
                'note':'Published mixed-model coefficients across pathology pseudo-progression in a Micro-PVM state; '
                       'not recomputed pseudobulk, not an identical Mathys contrast, not a rescue assay.'},
        protocol_hash=protocol['content_hash'], genes=sorted(genes,key=lambda g:g['gene_symbol']),
        counts={'up':counts[1],'down':counts[-1],'unmapped':len(unmapped)}, unmapped_symbols=unmapped,
        import_quality={'rows':len(table),'finite_p_family':int(finite.sum()),
                        'missing_p':int((~finite).sum()),'missing_effect':int(effects.isna().sum()),
                        'selected_before_mapping_and_cap':len(candidates)},
        l1000_note='External disease coefficients queried against the same binary LINCS context panel.'))
    # A sparse signature is a recorded outcome, never a reason to relax selection.
    if min(counts.values()) >= protocol['selection']['min_genes_per_direction']:
        validate_signature(bundle)
    return bundle


def compare_rankings(baseline, external):
    """Full union, exact compound IDs; missing coverage is never a zero score."""
    left = {r['pert_id']:r for r in baseline['ranking']}
    right = {r['pert_id']:r for r in external['ranking']}
    records = []; counts = {}
    for key in sorted(left.keys() | right.keys()):
        a, b = left.get(key), right.get(key)
        available = [r is not None and r.get('status') == 'scored' and
                     r.get('connectivity_mean') is not None for r in [a,b]]
        sa = a['connectivity_mean'] if available[0] else None
        sb = b['connectivity_mean'] if available[1] else None
        if not all(available): status = 'missing_coverage'
        elif sa == 0 or sb == 0: status = 'neutral'
        elif sa * sb < 0: status = 'opposite_direction'
        elif sa > 0: status = 'both_positive'
        else: status = 'both_negative'
        counts[status] = counts.get(status,0) + 1
        records.append({'pert_id':key, 'name':(a or b).get('pert_iname',key), 'comparison':status,
                        'baseline_score':sa,'external_score':sb,
                        'baseline_status':a['status'] if a else 'not_present',
                        'external_status':b['status'] if b else 'not_present',
                        'baseline_fdr':a.get('fdr') if available[0] else None,
                        'external_fdr':b.get('fdr') if available[1] else None})
    return {'counts':counts,'records':records,'interpretation':'Descriptive direction agreement, not a replication significance test.'}
