import numpy as np
import mdconf as md

file_name = 'conf_q5_walls.gro'

N = md.count_line(file_name)

input_file = open(file_name, 'r')

Z_half = 0.5*30.36600

n = 0
z_si_b = []
z_o1_b = []
z_o2_b = []
z_si_t = []
z_o1_t = []
z_o2_t = []

for line in input_file:
    if n>2 and n<N-1 :
        line_data = md.read_gro_line(line)
        if line_data[1] == "SUB" :
            if line_data[6] > Z_half :
                if line_data[2] == "SI":
                    z_si_t.append(line_data[6])
                elif line_data[2] == "O1":
                    z_o1_t.append(line_data[6])
                else:
                    z_o2_t.append(line_data[6])
            else :
                if line_data[2] == "SI":
                    z_si_b.append(line_data[6])
                elif line_data[2] == "O1":
                    z_o1_b.append(line_data[6])
                else:
                    z_o2_b.append(line_data[6])
    n += 1

input_file.close()

z_si_t = np.array(z_si_t)
z_o1_t = np.array(z_o1_t)
z_o2_t = np.array(z_o2_t)
z_si_b = np.array(z_si_b)
z_o1_b = np.array(z_o1_b)
z_o2_b = np.array(z_o2_b)

print("<z(O1)_bot> = "+str(np.mean(z_o1_b))+"nm")
print("<z(Si)_bot> = "+str(np.mean(z_si_b))+"nm")
print("<z(O2)_bot> = "+str(np.mean(z_o2_b))+"nm")
print("<z(O1)_top> = "+str(np.mean(z_o1_t))+"nm")
print("<z(Si)_top> = "+str(np.mean(z_si_t))+"nm")
print("<z(O2)_top> = "+str(np.mean(z_o2_t))+"nm")

print("Quad lenght top = "+str(np.abs(np.mean(z_o1_t)-np.mean(z_o2_t)))+"nm")
print("Quad lenght bot = "+str(np.abs(np.mean(z_o1_b)-np.mean(z_o2_b)))+"nm")
