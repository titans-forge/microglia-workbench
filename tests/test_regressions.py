"""Independent correctness cases and fail-closed output contracts."""
import itertools
import json

import numpy as np
import pandas as pd
import pytest

from microglia_repurposing.contexts import ContextKS, weighted_connectivity
from microglia_repurposing.ingest import parse_gmt
from microglia_repurposing.context_engine import _jsonable, _resolve_names
from microglia_repurposing.packet import build_packet
from microglia_repurposing.scoring import cosine_reversal, permutation_null, cmap_reversal
from microglia_repurposing.identity import pin_symbol
from microglia_repurposing.provenance import verify_source


def test_gmt_preserves_blank_description_and_rejects_duplicates(tmp_path):
    path = tmp_path / 'test.gmt'
    path.write_text('blank\t\tAPOE\tTREM2\tAPOE\nwithdesc\tdescription\tC1QA\tC1QB\n')
    assert parse_gmt(path) == {'blank': ['APOE', 'TREM2'], 'withdesc': ['C1QA', 'C1QB']}
    path.write_text('same\t\tA\nsame\t\tB\n')
    with pytest.raises(ValueError, match='duplicate_term'):
        parse_gmt(path)


def test_sparse_reversal_mimic_and_empty_have_known_answers():
    d = pd.Series([3, 2, 1, -1, -2, -3], index=list('ABCDEF'))
    vec = np.array([-1, -1, 0, 0, 0, 1])
    ctx = pd.DataFrame({'reverse': vec, 'mimic': -vec, 'empty': vec * 0, 'all': np.ones(6)}, index=d.index)
    assert ContextKS(d, ctx).scores(d) == pytest.approx([1, -1, 0, 0])


def test_vector_matches_independent_dense_walk_exhaustively_with_ties():
    # All 729 sign patterns, including no-hit/all-hit/one-sided signatures.
    ctx = pd.DataFrame(np.array(list(itertools.product([-1, 0, 1], repeat=6))).T, index=list('ABCDEF'))
    for weights in ([3, 2, 1, -1, -2, -3], [2, 2, 1, -1, -2, -2], [1]*6):
        d = pd.Series(weights, index=ctx.index)
        actual = ContextKS(d, ctx).scores(d)
        assert np.allclose(actual, weighted_connectivity(d, ctx), atol=1e-12)
        assert np.allclose(actual, ContextKS(d.iloc[::-1], ctx.iloc[::-1]).scores(d.iloc[::-1]))
        assert np.allclose(actual, -ContextKS(d, -ctx).scores(d))


def test_nonfinite_inputs_and_invalid_permutations_rejected():
    d = pd.Series([1., -1.], index=['A', 'B'])
    r = pd.DataFrame([[np.inf, 1]], columns=d.index)
    with pytest.raises(ValueError, match='nonfinite'):
        cosine_reversal(d, r)
    with pytest.raises(ValueError, match='n_perm'):
        permutation_null(d, r, pd.Series([0.]), n_perm=0)
    finite = pd.DataFrame([[-1e300, 1e300], [0, 0]], columns=d.index)
    assert cosine_reversal(d*1e300, finite).tolist() == pytest.approx([1, 0])


def test_permutation_observed_order_is_preserved():
    d = pd.Series([1., -2., 3., -4.], index=list('ABCD'))
    r = pd.DataFrame([[-1, 1, -1, 1], [1, -1, 1, -1]], columns=d.index, index=['a','b'])
    obs = cosine_reversal(d, r)
    forward = permutation_null(d, r, obs, n_perm=50)
    reverse = permutation_null(d, r, obs.iloc[::-1], n_perm=50)
    assert forward.sort_index().equals(reverse.sort_index())


def test_context_packet_retains_actual_statistic_and_flags():
    disease = {'genes': [], 'content_hash': 'hash', 'cell_context': 'human_brain_microglia', 'cell_state': 'test'}
    ranked = pd.DataFrame([{'status':'scored', 'connectivity_mean':0.6, 'cosine_mean':0.4, 'fdr':0.1,
                           'pubchem_inchikey':'known', 'tox_flag':True, 'pert_iname':'example'}], index=['BRD-1'])
    d = pd.Series([1.,-1.],index=['A','B'])
    matrix = pd.DataFrame([[-1.,1.]],index=ranked.index,columns=d.index)
    packet = build_packet(disease, ranked, matrix, d, mode='context')
    nxt = packet['next_experiment']
    assert packet['drug_cell_context'] == 'lincs_l1000_chem_pert_instance'
    assert nxt['score']['value'] == 0.6 and nxt['cosine_reversal'] == 0.4
    assert nxt['cmap_reversal'] is None
    assert 'toxic_moa_flag' in nxt['blocking_reasons']
    assert _jsonable(ranked)[0]['pert_id'] == 'BRD-1'
    assert _jsonable(ranked)[0]['pert_iname'] == 'example'
    assert build_packet(disease, ranked.iloc[:0], matrix, d, mode='context')['next_experiment'] is None
    json.dumps(packet, allow_nan=False)


def test_names_are_case_insensitive():
    wanted, missing = _resolve_names(['Vorinostat'], {'BRD-1':'vorinostat'})
    assert wanted == {'BRD-1'} and missing == []


def test_ambiguous_alias_is_not_arbitrarily_assigned(tmp_path):
    path = tmp_path / 'hgnc.tsv'
    path.write_text('symbol\thgnc_id\tstatus\talias_symbol\nA\tHGNC:1\tApproved\tX\nB\tHGNC:2\tApproved\tX|A\n')
    assert pin_symbol('X', hgnc_path=path) is None
    assert pin_symbol('A', hgnc_path=path)['hgnc_id'] == 'HGNC:1'


def test_source_hash_mismatch_is_detected(tmp_path):
    path = tmp_path / 'source'
    path.write_text('changed')
    path.with_name('source.provenance.json').write_text(json.dumps({'sha256': '0'*64}))
    with pytest.raises(ValueError, match='source_hash_mismatch'):
        verify_source(path)
