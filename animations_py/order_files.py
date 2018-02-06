import os
os.chdir('new_an3/')
files = os.listdir()
print(files)
num_files = []

for file in files:
	num = file[3:]
	num = num.split('.')[0]
	num_files.append(int(num))

sort_idx = sorted(range(len(num_files)),key=lambda k: num_files[k])

new_files = []
for i in range(len(sort_idx)):
	new_files.append(files[sort_idx[i]])

files = new_files
for i in range(len(files)):
	new_name = 'pic' + str(i) + '.png'
	os.rename(files[i],new_name)