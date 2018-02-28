from aquaduct.traj.sandwich import FramesRangeCollection, SmartRangeIncrement

frc = FramesRangeCollection()

frc.append(SmartRangeIncrement(5,10))
frc.append(SmartRangeIncrement(2,10))
frc.append(SmartRangeIncrement(3,10))
frc.append(SmartRangeIncrement(20,10))
frc.append(SmartRangeIncrement(12,5))
frc.append(SmartRangeIncrement(40,10))
frc.append(SmartRangeIncrement(3,100))

