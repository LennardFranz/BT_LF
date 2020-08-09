
import iotool as io

class TLeCroy:

  def __init__(self, fname, debug=None, raw=None):
    self._FileName = fname
    if raw:
      self._Raw = raw
    else:
      self._Raw = False
    self._Info = dict() 
    self._Header = dict() 
    self._x, self._y, self._Header, self._Info = io.readTrc(self._FileName, raw = self._Raw)
    if debug:
      self._Debug = debug
    else:
      self._Debug = False
  

  def PrintPrivate(self):
    print(" *** INFO::lecroy:PrintPrivate ***")
    print("  self._FileName        " + str(self._FileName))
    print("  self._Debug           " + str(self._Debug))
    print("  self._Raw             " + str(self._Raw))
    print("  ---")
    print("   self._Header")
    for item in self._Header.items():
      print(item)
    print("  ---")
    print("  self._Info")
    for item in self._Info.items():
      print(item)
    print(" ------------------------")
    if self._Debug:
      print(type(self._y[0]))


  def GetTrace(self):
    return self._x, self._y


  def GetY(self):
    return self._y


  def GetHeader(self):
    return self._Header


  def GetInfo(self):
    return self._Info


  def PlotTrace(self, raw=True):
    import matplotlib.pyplot as plt
    if raw:
      plt.plot(self._x, self._y, 'k', label=self._FileName)
      plt.ylabel('amplitude (a.u.)')
      plt.xlabel('time (a.u.)')
    else:
      milli = 1/1000.0
      plt.plot(self._x, self._y / milli, 'k', label=self._FileName)
      plt.ylabel('volts (mV)')
      plt.xlabel('time (s)')

    plt.title(self._FileName)
    plt.show()


