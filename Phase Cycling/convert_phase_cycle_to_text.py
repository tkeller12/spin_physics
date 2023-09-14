import numpy as np


#out = np.loadtxt('6-pulse_DQC_phase_cycle.csv', delimiter = ',', usecols = [1,2,3,4,5,6,9])
out = np.loadtxt('6-pulse_DQC_phase_cycle2.csv', delimiter = ',')


ph_1 = out[:,1].squeeze().astype(int)
ph_2 = out[:,2].squeeze().astype(int)
ph_3 = out[:,3].squeeze().astype(int)
ph_4 = out[:,4].squeeze().astype(int)
ph_5 = out[:,5].squeeze().astype(int)
ph_6 = out[:,6].squeeze().astype(int)

rx = out[:,7].squeeze()


ph_1_str = 'ph1 = [' + ', '.join(map(str,ph_1)) + ']\n'
ph_2_str = 'ph2 = [' + ', '.join(map(str,ph_2)) + ']\n'
ph_3_str = 'ph3 = [' + ', '.join(map(str,ph_3)) + ']\n'
ph_4_str = 'ph4 = [' + ', '.join(map(str,ph_4)) + ']\n'
ph_5_str = 'ph5 = [' + ', '.join(map(str,ph_5)) + ']\n'
ph_6_str = 'ph6 = [' + ', '.join(map(str,ph_6)) + ']\n'


det1_str = []
det2_str = []

det_dict = {0:('a','b'), 90:('-b', 'a'), 180: ('-a','-b'), 270: ('b', '-a')}

for ix in range(len(rx)):
    print(rx[ix])
    det = det_dict[rx[ix]]
    det1_str.append(det[0])
    det2_str.append(det[1])




det1_str = 'det1 = [' + ', '.join(det1_str) + ']\n'
det2_str = 'det2 = [' + ', '.join(det2_str) + ']\n'

print(det1_str)
print(det2_str)


with open('6-pulse_DQC_phase_cycle_20230914.txt', 'w') as f:
    f.write(ph_1_str)
    f.write(ph_2_str)
    f.write(ph_3_str)
    f.write(ph_4_str)
    f.write(ph_5_str)
    f.write(ph_6_str)

    f.write(det1_str)
    f.write(det2_str)




