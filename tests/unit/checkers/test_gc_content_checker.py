import pytest
from genedesign.checkers.gc_content_checker import GCContentChecker

@pytest.fixture
def gc_checker():
    checker = GCContentChecker()
    checker.initiate()
    return checker

def test_passing_gc(gc_checker):
    """
    A sequence with ~50% GC content should pass globally and locally.
    """
    # 50% GC, long enough for local windows
    cds = "ATGC" * 30  # 120 bp, 50% GC
    passed, gc, msg = gc_checker.run(cds)
    assert passed == True
    assert abs(gc - 0.5) < 0.01
    assert msg == "ok"

def test_too_low_global_gc(gc_checker):
    """
    A sequence with very low GC content should fail the global check.
    """
    cds = "ATTT" * 30  # 25% GC
    passed, gc, msg = gc_checker.run(cds)
    assert passed == False
    assert gc < 0.40

def test_too_high_global_gc(gc_checker):
    """
    A sequence with very high GC content should fail the global check.
    """
    cds = "GCCC" * 30  # 75% GC
    passed, gc, msg = gc_checker.run(cds)
    assert passed == False
    assert gc > 0.65

def test_failing_local_window(gc_checker):
    """
    A sequence that passes globally but has a local window with extreme GC should fail.
    """
    # First 50 bp all AT (0% GC), rest balanced to bring global GC into range
    low_gc_region = "AT" * 25          # 50 bp, 0% GC
    balanced_region = "ATGC" * 50     # 200 bp, 50% GC
    cds = low_gc_region + balanced_region
    passed, gc, msg = gc_checker.run(cds)
    assert passed == False

def test_empty_sequence(gc_checker):
    """
    An empty sequence should pass without error.
    """
    passed, gc, msg = gc_checker.run("")
    assert passed == True
    assert gc == 0.0
