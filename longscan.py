import u6
import LJScope3
import matplotlib.pyplot as plt

DAQ = u6.U6()               #u6 object named daq
DAQ.getCalibrationData()
Signal = DAQ.getAIN(0)      #reads from channel 0

scanner = LJScope3.LJScan(DAQ)
scanner.set_rate(48_000)
scanner.set_n_scan(4000)       #scans n * 1200, so in this case 4800
times4,vals4 = scanner.scan()
scanner.close()

plt.plot(times4, vals4, 'r.')
