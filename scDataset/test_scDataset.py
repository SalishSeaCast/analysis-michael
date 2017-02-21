import os,fnmatch
from scDataset import scDataset
from IPython import embed

def gethourlyfiles(prefix, days, grid):
    files = []    
    for day in days:
        dirname = os.path.join(prefix, day)
        for item in os.listdir(dirname):
            if fnmatch.fnmatchcase(item, "SalishSea_1h_*"+grid+"*.nc"):
                files += [os.path.join(dirname,item)]
    return files

prefix = '/results/SalishSea/nowcast-blue'
days = ['{:02d}dec16'.format(x) for x in range(1,31+1)]
files = gethourlyfiles(prefix, days, 'grid_T')

print("Found {} files".format(len(files)))

print("Test 1")
ds = scDataset(files)
t = ds.variables['time_counter'][:]
ds.close()

print("Test 2")
ds = scDataset(files)
t1 = ds.variables['votemper'][:,:,100,100]
ds.close()

print("Test 3")
ds = scDataset(files)
for ti in range(0,3):
    print("Loading time "+str(ti))
    surfsal = ds.variables['vosaline'][ti,0,:,:]
    # make_a_plot(surfsal)
ds.close()

print("Test 4")
ds = scDataset(files)
t2 = ds.variables['votemper'][5:-1,1:8,100:190,100:-1]
print(t2.shape)
ds.close()

#embed()
