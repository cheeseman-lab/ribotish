"""Process pools for Ribo-TISH (fork): non-daemonic outer pool, ordered results.

Added 2026-10-06 by Matteo DiBernardo for the Ribo-TISH Cheeseman lab fork. Part of Ribo-TISH
(Copyright (c) 2016 Peng Zhang), distributed under the GNU General Public License v3 (LICENSE).
"""
import multiprocessing
import multiprocessing.pool

_CTX = multiprocessing.get_context("fork")


class NoDaemonProcess(_CTX.Process):
  @property
  def daemon(self):
    return False

  @daemon.setter
  def daemon(self, value):
    pass


class NoDaemonContext(type(_CTX)):
  Process = NoDaemonProcess


def nodaemon_pool(processes):
  """Pool whose workers may start their own pools (used for TIS background estimation)."""
  return multiprocessing.pool.Pool(processes, context=NoDaemonContext())


def pool(processes):
  return _CTX.Pool(processes)


def _call_trapping_exit(func, args, kwargs):
  # SystemExit kills a pool worker without reporting back, so apply() would wait forever.
  try:
    return False, func(*args, **kwargs)
  except SystemExit as e:
    return True, e.code


def run_in_child(func, *args, **kwargs):
  """Run func in one non-daemonic child process (frees its memory on return); re-raise exit()."""
  p = nodaemon_pool(1)
  try:
    exited, value = p.apply(_call_trapping_exit, (func, args, kwargs))
  finally:
    p.close()
    p.join()
  if exited:
    raise SystemExit(value)
  return value
