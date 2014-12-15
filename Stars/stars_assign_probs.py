import numpy as np

name_dict= {'Andrew E.':[], 'Sarah':[], 'Andrea':[], 'Adrian':[],\
            'Jingjing': [], 'David':[], 'Emily':[], 'Jose':[],
            'Andrew W.':[]}

#\
 #        'Jose':[]  , 'Andrew W.':[], 'Emily': []




prob_nums = range(1,101,1)

# doing some silly surgery
surgery = [86]

print np.where(np.array(prob_nums) == 86)[0][0]

for s in surgery:
#    print prob_nums[np.where(prob_nums) == s)]
 #   print prob_nums[s]
    prob_nums.pop( np.where(np.array(prob_nums)==s)[0][0] )


l = len(prob_nums)



i = 0

keys = name_dict.keys()
key_len = len(keys)
while l > 0:
  
    name = keys[i]
    
    index = np.random.randint(0, l) # pick an index

    name_dict[name].append(prob_nums[index]) # assign prob number
    prob_nums.pop(index)                     # remove prob number

    l = len(prob_nums)
    i = i + 1           # written so one number per person at a time
    if (i == key_len): 
        i = 0           # time to start over         
        np.random.shuffle(keys) # randomly shuffle order of people

    
for name in name_dict:
    
    name_dict[name].sort() # sort from low to high 

    print "%9s"%(name), "[ ",
    len_array = len(name_dict[name])
    for i in range(len_array):
    
        print "%3i"%(name_dict[name][i]),

    
    print " ]", len_array
