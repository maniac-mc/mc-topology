# water "SPC/E" parameters ------------------------------
bond_coeff   1     600.0   1.0 
angle_coeff  1   75.0    109.47
pair_coeff   1 1  0.1553  3.166 
pair_coeff   2 2  0.0     0.0
group O type 1
group H type 2
group water type  1  2
fix fShakeSPCE water shake 0.0001 20 0 b 1 a 1  
# (Remember to "unfix" fShakeSPCE during minimization.)
# Lennard Jones wall parameters -------------------------
pair_coeff  3 3 0.164     3.374
group wall type 3
