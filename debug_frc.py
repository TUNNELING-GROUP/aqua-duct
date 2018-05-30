from aquaduct.traj.sandwich import FramesRangeCollection
from aquaduct.utils.helpers import lind, SmartRangeIncrement

frc = FramesRangeCollection()

frc.append(SmartRangeIncrement(5,10))
frc.append(SmartRangeIncrement(2,10))
frc.append(SmartRangeIncrement(3,10))
frc.append(SmartRangeIncrement(20,10))
frc.append(SmartRangeIncrement(12,5))
frc.append(SmartRangeIncrement(40,10))
frc.append(SmartRangeIncrement(3,100))


l = list(frc.get_ranges(SmartRangeIncrement(12,30)))

s = list(SmartRangeIncrement(12,30).get())

ss = []
for sr,xr in l:
    ss.append(lind(list(sr.get()),xr))




