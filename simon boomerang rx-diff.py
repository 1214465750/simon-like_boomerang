import math

import numpy as np


'''
simon  block    key
         32     64
         48     72 96
         64    96  128
         96     96 144
         128  128 192 256
'''
'''
here to choose the version of simeck  : 
list_version[0] is simon32/64  
list_version[1] is simon48/72
list_version[2] is simon48/96
list_version[3] is simon64/96
list_version[4] is simon64/128
list_version[5] is simon96/96
list_version[6] is simon96/144
list_version[7] is simon128/128  
list_version[8] is simon128/192  
list_version[9] is simon128/256      
'''
list_version = [[16, 4, 32, 0], [24, 3, 36, 0], [24, 4, 36, 1], [32, 3, 42, 2], [32, 4, 44, 3], [48, 2, 52, 2],
                    [48, 3, 54, 3], [64, 2, 68, 2], [64, 3, 69, 3], [64, 4, 72, 4]]
block_n, num_of_key, total_round, z_num = list_version[0]

date_num = 2 ** 15
date_block=2**12
block_size = 2 * block_n
num_of_key = 4
key_size = num_of_key * block_n
'''
round: the round of distinguisher
diff_in_l_int: left part of input difference of E0
diff_in_r_int: right part of input difference of E0
diff_out_l_int: left part of output difference of E1
diff_out_r_int: right part of output difference of E1
diff_mainkey_1: mainkey difference of E0
diff_mainkey_2: last four round key difference of E1
'''
round = 14

diff_in_l_int=0x6
diff_in_r_int=0x617
diff_out_l_int=0x8000
diff_out_r_int=0x0
diff_mainkey_1=[0x3 ,0x6 ,0x0 ,0x0 ,]
diff_mainkey_2=[0x1f58 ,0x8880 ,0xe540 ,0x880 ,]
offset=1

z = np.zeros((5, 62), dtype=bool)
const_c = np.zeros(block_n, dtype=bool)
for i in range(2, block_n):
    const_c[i] = 1
z = [[1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
      0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, ],
     [1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1,
      0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, ],
     [1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, ],
     [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
      0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, ],
     [1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, ]]


def sl_1(a, b):

    b=b%block_n
    c=a.copy()
    c[:b]=a[block_n-b:]
    c[b:] = a[:block_n-b]
    return c

def sl_2(a, b):

    b=b%block_n
    c=a.copy()
    c[:,:b]=a[:,block_n-b:]
    c[:, b:] = a[:,:block_n-b]
    return c

def sl_3(a, b):

    b=b%block_n
    c=a.copy()
    c[:,:,:b]=a[:,:,block_n-b:]
    c[:,:, b:] = a[:,:,:block_n-b]
    return c


def keyschedule(mainkey, round):
    k = np.zeros((date_block,round, block_n), dtype=bool)
    k[:,:num_of_key] = mainkey

    for i in range(num_of_key, round):
        tmp = sl_2(k[:,i - 1], -3)
        if num_of_key == 4:
            tmp ^= k[:,i - 3]
        tmp = tmp ^ sl_2(tmp, -1)
        k[:,i] = k[:,i - num_of_key] ^ tmp ^ const_c
        k[:,i,0] ^= bool(z[z_num][(i - num_of_key) % 62])



    return k


def ni_keyschedule(mainkey,round):
    k = np.zeros((date_block,round, block_n), dtype=bool)
    k[:,round-num_of_key:round] = mainkey

    for i in range(round-num_of_key-1,-1,-1):
        tmp = sl_2(k[:,i + num_of_key-1], -3)
        if num_of_key == 4:
            tmp ^= k[:,i +1]
        tmp = tmp ^ sl_2(tmp, -1)
        k[:,i] = k[:,i + num_of_key] ^ tmp ^ const_c
        k[:,i, 0] ^= bool(z[z_num][(i) % 62])

    return k

def zhong_keyschedule(mainkey,round,zhong):
    k = np.zeros((date_block,round, block_n), dtype=bool)
    k[:,zhong:zhong+num_of_key]=mainkey

    for i in range(zhong-1,-1,-1):
        tmp = sl_2(k[:,i + num_of_key - 1], -3)
        if num_of_key == 4:
            tmp ^= k[:,i + 1]
        tmp = tmp ^ sl_2(tmp, -1)
        k[:,i] = k[:,i + num_of_key] ^ tmp ^ const_c
        k[:,i, 0] ^= bool(z[z_num][(i) % 62])


    for i in range(zhong+num_of_key,round):
        tmp = sl_2(k[:,i - 1], -3)
        if num_of_key == 4:
            tmp ^= k[:,i - 3]
        tmp = tmp ^ sl_2(tmp, -1)
        k[:,i] = k[:,i - num_of_key] ^ tmp ^ const_c
        k[:,i, 0] ^= bool(z[z_num][(i - num_of_key) % 62])

    return k
def oneenc(plain_l, plain_r, key):
    tmp = plain_l
    plain_l = plain_r ^ (sl_2(plain_l, 1) & sl_2(plain_l, 8)) ^ sl_2(plain_l, 2) ^ key
    plain_r = tmp
    return plain_l, plain_r


def enc(plain_l, plain_r, roundkey, round):

    for i in range(round):
        plain_l, plain_r = oneenc(plain_l, plain_r, roundkey[:,i])
    return plain_l, plain_r


def onedec(plain_l, plain_r, key):

    tmp = plain_r
    plain_r = plain_l ^ (sl_2(plain_r, 1) & sl_2(plain_r, 8)) ^ sl_2(plain_r, 2) ^ key
    plain_l = tmp
    return plain_l, plain_r





def dec(ciper_l, ciper_r, roundkey, round):
    for i in range(round):
        ciper_l, ciper_r = onedec(ciper_l, ciper_r, roundkey[:,round-1-i])

    return ciper_l, ciper_r


def array_to_int(a):
    b = 0
    for i in range(len(a)):
        b += a[i] << i
    return b

def int_to_array(a):
    arr=np.zeros(block_n, dtype=bool)
    for i in range(block_n):
        arr[i]=(a>>i)%2
    return arr
def veirfy_1(round, diff_in_l, diff_in_r, diff_out_l, diff_out_r,diff_mainkey_1,diff_mainkey_2):



    sum = 0
    for i in range(date_num//date_block):
        mainkey_1 = np.random.randint(0, 2, size=(date_block,num_of_key, block_n))
        t = sl_3(mainkey_1,offset)
        mainkey_2 = t[:]^diff_mainkey_1
        mainkey_3 = mainkey_1[:] ^ diff_mainkey_2
        mainkey_4 = mainkey_2[:] ^ sl_2(diff_mainkey_2,offset)
        roundkey_1 = keyschedule(mainkey_1, round)
        roundkey_2 = keyschedule(mainkey_2, round)
        roundkey_3 = keyschedule(mainkey_3, round)
        roundkey_4 = keyschedule(mainkey_4, round)

        plain_l_1 = np.random.randint(0, 2, size=(date_block,block_n))
        plain_r_1 = np.random.randint(0, 2, size=(date_block,block_n))
        plain_l_2 = sl_2(plain_l_1,offset) ^ diff_in_l
        plain_r_2 = sl_2(plain_r_1,offset) ^ diff_in_r
        ciper_l_1, ciper_r_1 = enc(plain_l_1, plain_r_1, roundkey_1, round)
        ciper_l_2, ciper_r_2 = enc(plain_l_2, plain_r_2, roundkey_2, round)
        ciper_l_3 = ciper_l_1 ^ diff_out_l
        ciper_r_3 = ciper_r_1 ^ diff_out_r
        ciper_l_4 = ciper_l_2 ^ sl_1(diff_out_l,offset)
        ciper_r_4 = ciper_r_2 ^ sl_1(diff_out_r,offset)
        plain_l_3, plain_r_3 = dec(ciper_l_3,ciper_r_3,roundkey_3, round)
        plain_l_4, plain_r_4 = dec(ciper_l_4, ciper_r_4, roundkey_4, round)

        sum+=np.sum(~((sl_2(plain_l_3,offset) ^ plain_l_4^diff_in_l).any(axis=1)|(sl_2(plain_r_3,offset) ^ plain_r_4 ^ diff_in_r).any(axis=1)))
    print(sum)
    print(date_num)
    if sum!=0:
        print(math.log2(date_num / sum))

    return None


def sp(a):
    print(np.shape(a))

if __name__ == '__main__':




    diff_in_l_arr = int_to_array(diff_in_l_int)
    diff_in_r_arr = int_to_array(diff_in_r_int)
    diff_out_l_arr = int_to_array(diff_out_l_int)
    diff_out_r_arr = int_to_array(diff_out_r_int)
    diff_mainkey_1_arr=np.array([int_to_array(diff_mainkey_1[i])for i in range(4)])
    diff_mainkey_2_arr = np.array([int_to_array(diff_mainkey_2[i]) for i in range(4)])

    veirfy_1(round, diff_in_l_arr, diff_in_r_arr, diff_out_l_arr, diff_out_r_arr,diff_mainkey_1_arr,diff_mainkey_2_arr)





