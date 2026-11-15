"""Tests for the multiple sequence aligner selection (FAMSA default, MUSCLE optional)."""
import os
import shutil
import stat

import pytest

from bofasa import alignment, config
from bofasa.main import create_main_parser


FASTA = ">s1|a\nMKVLAAGIVGLLLA\n>s2|b\nMKVLAGGIVGLLA\n>s3|c\nMKVIAAGIVGLLLA\n"


def _write(path, text):
    path.write_text(text)
    return str(path)


def _fake_aligner(bin_dir, name, body):
    """Create an executable stand-in for an aligner on PATH."""
    exe = bin_dir / name
    exe.write_text("#!/bin/sh\n" + body)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return exe


def test_default_aligner_is_famsa():
    assert config.DEFAULT_ALIGNER == "famsa"
    assert config.ALIGNMENT_METHODS == ["famsa", "muscle"]


def test_build_alignment_cmd():
    assert alignment.build_alignment_cmd("in.faa", "out.faa", 4) == [
        "famsa", "-t", "4", "in.faa", "out.faa"]
    assert alignment.build_alignment_cmd("in.faa", "out.faa", 2, "muscle") == [
        "muscle", "-super5", "in.faa", "-threads", "2", "-output", "out.faa"]
    with pytest.raises(ValueError):
        alignment.build_alignment_cmd("in.faa", "out.faa", 1, "clustalo")


@pytest.mark.parametrize("aligner", ["famsa", "muscle"])
def test_run_alignment_uses_requested_aligner(tmp_path, monkeypatch, aligner):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "called.txt"
    # famsa: <in> is $3, <out> is $4; muscle: <in> is $2, <out> is $6.
    io_args = '"$3" "$4"' if aligner == "famsa" else '"$2" "$6"'
    _fake_aligner(bin_dir, aligner,
                  f'echo "$0 $@" > "{log}"\nset -- {io_args}\ncp "$1" "$2"\n')
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")

    in_faa = _write(tmp_path / "og.faa", FASTA)
    out_faa = str(tmp_path / "og.msa.faa")
    assert alignment.run_alignment([in_faa, out_faa, 3, aligner]) is True
    called = log.read_text()
    assert os.path.basename(called.split()[0]) == aligner
    assert os.path.isfile(out_faa)


def test_run_alignment_defaults_to_famsa_for_three_item_jobs(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _fake_aligner(bin_dir, "famsa", 'cp "$3" "$4"\n')
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    in_faa = _write(tmp_path / "og.faa", FASTA)
    assert alignment.run_alignment([in_faa, str(tmp_path / "out.faa"), 1]) is True


def test_run_alignment_reports_failure(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _fake_aligner(bin_dir, "famsa", "exit 3\n")
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    in_faa = _write(tmp_path / "og.faa", FASTA)
    assert alignment.run_alignment([in_faa, str(tmp_path / "out.faa"), 1, "famsa"]) is False


def test_run_alignment_skips_single_sequence(tmp_path):
    in_faa = _write(tmp_path / "og.faa", ">s1|a\nMKV\n")
    out_faa = tmp_path / "out.faa"
    assert alignment.run_alignment([str(in_faa), str(out_faa), 1, "famsa"]) is True
    assert not out_faa.exists()


@pytest.mark.skipif(shutil.which("famsa") is None, reason="famsa not installed")
def test_run_alignment_with_real_famsa(tmp_path):
    in_faa = _write(tmp_path / "og.faa", FASTA)
    out_faa = str(tmp_path / "og.msa.faa")
    assert alignment.run_alignment([in_faa, out_faa, 1, "famsa"]) is True


@pytest.mark.skipif(shutil.which("muscle") is None, reason="muscle not installed")
def test_run_alignment_with_real_muscle(tmp_path):
    in_faa = _write(tmp_path / "og.faa", FASTA)
    out_faa = str(tmp_path / "og.msa.faa")
    assert alignment.run_alignment([in_faa, out_faa, 1, "muscle"]) is True


def test_cli_aligner_option():
    parser = create_main_parser()
    args = parser.parse_args(["run", "-i", "prep", "-o", "out"])
    assert args.aligner == "famsa"
    args = parser.parse_args(["run", "-i", "prep", "-o", "out", "--aligner", "muscle"])
    assert args.aligner == "muscle"
    args = parser.parse_args(["run", "-i", "prep", "-o", "out", "-al", "famsa"])
    assert args.aligner == "famsa"
    with pytest.raises(SystemExit):
        parser.parse_args(["run", "-i", "prep", "-o", "out", "--aligner", "clustalo"])
