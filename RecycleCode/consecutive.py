import numpy as np

nums = np.array([1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,1,0,0,0,1]).astype(bool)

def getIndexMaxOnes(arr):
    ## returns the starting index and ending index of the max consecutive ones
    # intitialize count
    cur_count = 0
    cur_sta = 0

    max_count = 0
    pre_state = 0

    index_sta = 0
    index_end = 0

    for i in range(0, np.size(arr)):

        if (arr[i] == 0):
            cur_count = 0
            if((pre_state == 1)&(cur_sta == index_sta)):
                index_end = i-1
            pre_state = 0


        else:
            if(pre_state == 0):
                cur_sta = i
                pre_state = 1
            cur_count+= 1
            if(cur_count>max_count):
                max_count = cur_count
                index_sta = cur_sta

    return index_sta,index_end

index_sta, index_end = getIndexMaxOnes(nums)

print(index_sta, index_end)
