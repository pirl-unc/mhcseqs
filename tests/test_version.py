def test_version_importable():
    from mhcseqs.version import __version__

    assert isinstance(__version__, str)
    assert "." in __version__


def test_version_from_init():
    import mhcseqs

    assert hasattr(mhcseqs, "__version__")
    assert mhcseqs.__version__ == mhcseqs.version.__version__
