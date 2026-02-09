def test_import() -> None:
    try:
        import torsion_gertsenshtein

        assert torsion_gertsenshtein is not None, (
            "torsion_gertsenshtein module should not be None"
        )
    except ImportError as e:
        msg = "Failed to import torsion_gertsenshtein module"
        raise AssertionError(msg) from e
