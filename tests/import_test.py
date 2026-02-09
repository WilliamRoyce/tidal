def test_import() -> None:
    try:
        import tidal

        assert tidal is not None, (
            "tidal module should not be None"
        )
    except ImportError as e:
        msg = "Failed to import tidal module"
        raise AssertionError(msg) from e
