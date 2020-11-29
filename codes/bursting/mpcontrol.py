##provide ^c for multiprocessing

import multiprocessing
import psutil
import signal



class _multiprocessing_worker_init(object):
  def __init__(self, parent_id):
    self.parent_id = parent_id
  def __call__(self):
    def sig_int(signal_num, frame):
      print('signal: %s' % signal_num)
      parent = psutil.Process(self.parent_id)
      for child in parent.children():
        if child.pid != os.getpid():
          print("killing child: %s" % child.pid)
          child.kill()
      print("killing parent: %s" % self.parent_id)
      parent.kill()
      print("suicide: %s" % os.getpid())
      psutil.Process(os.getpid()).kill()
    signal.signal(signal.SIGINT, sig_int)