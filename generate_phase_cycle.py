

pulse = ['90','180']
cycle = [(0,90,180,270), (0,180)]
steps = [4,2]
delta_p = [1,2]
#rx_phase = [(a,b), (a,b)]

total = 1
for each in cycle:
    total *= len(each)

print('total', total)


complete = []

for 



for ix, name in enumerate(pulse):
    print(ix, name)

    for num in range(steps[ix]):
        print(num)


#def gen_phase_cycle(pulse, steps, delta_p):
#    if len(pulse) == 1:
#        return
#    else
#        gen_phase_cycle(pulse[1:],steps[1:],delta_p[1:])
#
#gen_phase_cycle(pulse, steps, delta_p)
