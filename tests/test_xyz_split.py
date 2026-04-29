"""Tests for the multi-frame XYZ splitter."""

from __future__ import annotations

from rdkit import Chem

from app.modules.toolkits.helpers import split_xyz_frames


WATER = (
    "3\n"
    "water\n"
    "O    0.0000    0.0000    0.0000\n"
    "H    0.7572    0.5860    0.0000\n"
    "H   -0.7572    0.5860    0.0000\n"
)
METHANE = (
    "5\n"
    "methane\n"
    "C    0.0000    0.0000    0.0000\n"
    "H    0.6300    0.6300    0.6300\n"
    "H   -0.6300   -0.6300    0.6300\n"
    "H   -0.6300    0.6300   -0.6300\n"
    "H    0.6300   -0.6300   -0.6300\n"
)


class TestSplitXYZFrames:
    def test_mid_stream_truncated_frame_is_dropped_and_next_frame_recovers(self):
        """Truncated frame mid-stream must not consume lines from the next frame."""
        truncated = "5\nbad\nC 0 0 0\nH 1 1 1\n"  # declares 5 atoms, only 2 supplied
        frames = split_xyz_frames(WATER + truncated + METHANE)
        assert len(frames) == 2
        comments = {f.splitlines()[1].strip() for f in frames}
        assert comments == {"water", "methane"}

    def test_atom_row_validation_rejects_consumed_next_frame_header(self):
        """A frame that declares 5 atoms but only has 2 real atom rows must be rejected.

        Without per-atom-row validation the splitter would greedily treat the
        next frame's count line and comment as pseudo-atom rows, packaging a
        bogus 5-row block instead of recovering the following valid frame.
        """
        # Frame A: declares 5 atoms, provides only 2 real atom rows; then Frame B follows.
        frame_a = "5\nbogus\nC 0 0 0\nH 1 1 1\n"
        frames = split_xyz_frames(frame_a + METHANE)
        # frame_a must be rejected; METHANE must be recovered.
        assert len(frames) == 1
        assert frames[0].splitlines()[1].strip() == "methane"

    def test_unicode_superscript_digit_not_treated_as_atom_count(self):
        """A line consisting of a Unicode superscript digit must not be parsed as a count."""
        superscript_line = "²\n"  # U+00B2 — looks like '2' but is not ASCII
        frames = split_xyz_frames(WATER + superscript_line + METHANE)
        assert len(frames) == 2
        comments = {f.splitlines()[1].strip() for f in frames}
        assert comments == {"water", "methane"}

    def test_padded_count_line_is_normalized(self):
        """The count line in the emitted frame must be stripped of surrounding whitespace."""
        padded = "  3  \nwater\nO 0 0 0\nH 1 0 0\nH -1 0 0\n"
        frames = split_xyz_frames(padded)
        assert len(frames) == 1
        assert frames[0].splitlines()[0] == "3"

    def test_empty_returns_empty_list(self):
        assert split_xyz_frames("") == []

    def test_whitespace_only_returns_empty_list(self):
        assert split_xyz_frames("   \n  \n") == []

    def test_single_frame(self):
        frames = split_xyz_frames(WATER)
        assert len(frames) == 1
        assert frames[0].splitlines()[0].strip() == "3"

    def test_two_frames_concatenated(self):
        frames = split_xyz_frames(WATER + METHANE)
        assert len(frames) == 2
        assert frames[0].splitlines()[1].strip() == "water"
        assert frames[1].splitlines()[1].strip() == "methane"

    def test_three_frames(self):
        frames = split_xyz_frames(WATER + METHANE + WATER)
        assert len(frames) == 3
        assert frames[0].splitlines()[1].strip() == "water"
        assert frames[1].splitlines()[1].strip() == "methane"
        assert frames[2].splitlines()[1].strip() == "water"

    def test_crlf_line_endings(self):
        frames = split_xyz_frames(WATER.replace("\n", "\r\n"))
        assert len(frames) == 1
        # Splitter normalizes to Unix line endings inside each frame.
        assert "\r\n" not in frames[0]

    def test_no_trailing_newline(self):
        frames = split_xyz_frames(WATER.rstrip())
        assert len(frames) == 1
        # Each emitted frame must end with exactly one newline.
        assert frames[0].endswith("\n")

    def test_truncated_last_frame_is_dropped(self):
        # Declares 5 atoms but supplies only 2.
        truncated = "5\nbad\nC 0 0 0\nH 1 1 1\n"
        frames = split_xyz_frames(WATER + truncated)
        assert len(frames) == 1  # only the WATER frame survives
        assert frames[0].splitlines()[1].strip() == "water"

    def test_garbage_between_frames_is_skipped(self):
        garbage = "garbage line 1\ngarbage line 2\n"
        frames = split_xyz_frames(WATER + garbage + METHANE)
        # Splitter advances past non-numeric lines and recovers the next frame.
        assert len(frames) == 2
        assert {f.splitlines()[1].strip() for f in frames} == {"water", "methane"}

    def test_unicode_in_comment_is_preserved(self):
        title = "ferrocène (η5)"
        frame = f"3\n{title}\nO 0 0 0\nH 0 0 1\nH 0 1 0\n"
        frames = split_xyz_frames(frame)
        assert frames[0].splitlines()[1] == title

    def test_atom_count_zero_yields_empty_frame_dropped(self):
        # Zero atoms is a degenerate XYZ; we drop it rather than emit a no-atom mol.
        frame = "0\nempty\n"
        frames = split_xyz_frames(frame + WATER)
        assert len(frames) == 1
        assert frames[0].splitlines()[1].strip() == "water"

    def test_returns_self_contained_frames(self):
        """Each returned string must itself be a valid XYZ block (parseable on its own)."""
        frames = split_xyz_frames(WATER + METHANE)
        for f in frames:
            mol = Chem.MolFromXYZBlock(f)
            assert mol is not None
            assert mol.GetNumAtoms() in (3, 5)
