def test_import() -> None:
    try:
        import torsion_gertsenshtein  # noqa: PLC0415
    except ImportError:
        torsion_gertsenshtein = None

    assert torsion_gertsenshtein is not None, "torsion_gertsenshtein module should not be None"
