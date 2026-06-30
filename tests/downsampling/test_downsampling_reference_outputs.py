from utils import (
    TEST_SAMPLES,
    assert_stability_table,
    reference_output_for_sample,
    repo_root,
)


def test_reference_outputs_exist_and_have_expected_schema():
    root = repo_root()
    for sample in TEST_SAMPLES:
        reference = root / reference_output_for_sample(sample)
        assert reference.exists(), f"missing reference output: {reference}"
        assert_stability_table(reference)
