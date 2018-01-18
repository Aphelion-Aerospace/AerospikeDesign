import numpy as np 

def clean_ID_list(initial_ID_list,del_ID_list):
    ID_list = initial_ID_list.copy()

    for ID in ID_list:
        if ID in del_ID:
            ID_list.remove(ID)

    ID_array = np.asarray(ID_list)
    ID_array_original = ID_array.copy()
    for ID in del_ID_list:
        ID_array[ID_array_original > ID] = ID_array[ID_array_original > ID] - 1
        print(ID_array)
    
    ID_list = []
    for ID in ID_array:
        ID_list.append(ID)

    return ID_list 

a = np.array([1,2,3,4,5,6,7,8,9,10]) # ie. a[2] = 3 etc...

b = [1,4,5,7,8,9] # ie. reference array(2,5,6,8,9,10)

del_ID = [3,6,7,9] # ie. delete 4,7,8,10

# after delete... 
#	a = [1,2,3,5,6,9]
#		[0,1,2,3,4,5]
#	a[b] = array(2,5,6,9),
#	b = [1,3,4,5]

print('Initially, a[b] = ' + str(a[b]))

a = np.delete(a,del_ID)

print("a = " + str(a) + '. Should be (1,2,3,5,6,9)')

b = clean_ID_list(b,del_ID)

print(b)

print("a[b] = " + str(a[b]) + '. Should be (2,5,6,9)')