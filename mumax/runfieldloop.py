import subprocess
import shutil
from tqdm import tqdm 

seed_start = 0
seed_end = 20

for i in tqdm(range(seed_start,seed_end)):
    file = open('fieldloop.mx3','r')
    content = file.readlines()
    content[18] = f'ThermSeed({i})\n'
    file = open('fieldloop.mx3','w')
    file.writelines(content)
    file.close()
    subprocess.run(['mumax3','fieldloop.mx3'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    shutil.copyfile(f'./fieldloop.out/table.txt',f'./fieldloop_r{i}.txt')