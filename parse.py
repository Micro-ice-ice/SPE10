layer = int(input("Enter layer "))
output = []
with open('spe_perm.dat') as file1:
    with open('spe_phi.dat') as file2:
        data1 = file1.read().splitlines()
        data2 = file2.read().splitlines()
        for line in data1:
            row = line.split('\t')
            for elem in row:
                value = elem.strip()
                if value == '':
                    continue
                output.append(value)
        for line in data2:
            row = line.split('\t')
            for elem in row:
                value = elem.strip()
                if value == '':
                    continue
                output.append(value)
file2.close()
file1.close()

kx = []
ky = []
phi = []
step = 60 * 220 * 85

for j in range(60 * 220 * layer, 60 * 220 * (layer + 1)):
    kx.append(output[j])
    ky.append(output[j + step])
    phi.append(output[j + 3 * step])

file = open("data.txt", "w")
for i in range(len(kx)):
    file.write(kx[i] + " " + ky[i] + " " + phi[i] + " ")
file.close()


