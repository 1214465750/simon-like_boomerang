import math
from time import time
import numpy as np


'''
simeck的差分链搜索  block    key
         32     64
         48     96
         64    128
'''
'''
round: the round of distinguisher
diff_in_l_int: left part of input difference of E0
diff_in_r_int: right part of input difference of E0
diff_out_l_int: left part of output difference of E1
diff_out_r_int: right part of output difference of E1
diff_mainkey_1: mainkey difference of E0
diff_mainkey_2: last four round key difference of E1
'''
round=14
diff_in_l_int=0x4
diff_in_r_int=0xc
diff_out_l_int=0xa
diff_out_r_int=0x4
diff_mainkey_1=[0x6 ,0x0 ,0x2 ,0x0 ,]
diff_mainkey_2=[0x0 ,0x2 ,0x0 ,0x0 ,]
'''
offset: rotation-xor offset 
'''
offset=1
date_num = 2 ** 22
date_block=2**19
'''
here to choose the version of simeck  : 
version[0] is simeck32/64  
version[1] is simeck48/96
version[2] is simeck64/128
'''
version=[[16,32,0],[24,36,0],[32,44,1],]
block_n,total_round,z_choose=version[0]
block_size = 2 * block_n
num_of_key = 4
key_size = num_of_key * block_n

z0 = [1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0]
z1 = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1,
      1,
      0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0]
z0 = np.array(z0, dtype=bool)
z1 = np.array(z1, dtype=bool)
z = [z0, z1]
a = 5
b = 0
c = 1
def sl(a, b):
    return np.roll(a, b)
def sl_3(a, b):
    b=b%block_n
    c=a.copy()
    c[:,:,:b]=a[:,:,block_n-b:]
    c[:,:, b:] = a[:,:,:block_n-b]
    return c


def sl_2(a, b):
    b=b%block_n
    c=a.copy()
    c[:,:b]=a[:,block_n-b:]
    c[:, b:] = a[:,:block_n-b]
    return c



def keyschedule(mainkey, round,block_n,num_of_key,a,b,c,z):
    k = np.zeros((date_block,round,block_n), dtype=bool)
    k[:,:num_of_key] = mainkey
    rc=np.ones((date_block,block_n),dtype=bool)
    rc[:,0]=0
    rc[:,1] = 0
    for i in range(num_of_key,round):
        rk=rc.copy()
        rk[:,0]=z[(i-4)%len(z)]
        k[:,i]=oneenc(k[:,i-3],k[:,i-4],rk,a,b,c,)
    return k


def ni_keyschedule(mainkey,round):
    k = np.zeros((date_block,round, block_n), dtype=bool)
    k[:,round-num_of_key:] = mainkey

    rc = np.ones((date_block, block_n), dtype=bool)
    rc[:, 0] = 0
    rc[:, 1] = 0
    for i in range(round-num_of_key-1,-1,-1 ):
        rk = rc.copy()
        rk[:, 0] = z[z_choose][i % len(z[z_choose])]
        k[:, i ] = onedec(k[:, i +4 ], k[:, i +1], rk, a, b, c, )
    return k

def oneenc(plain_l, plain_r, key,a,b,c):

    plain_l = plain_r ^ (sl_2(plain_l, a) & sl_2(plain_l, b)) ^ sl_2(plain_l, c) ^ key

    return plain_l


def enc(plain_l, plain_r, roundkey, round,a,b,c):

    for i in range(round):
        tmp=plain_l
        plain_l = oneenc(plain_l, plain_r, roundkey[:,i],a,b,c)
        plain_r=tmp
    return plain_l, plain_r


def onedec(plain_l, plain_r, key,a,b,c):


    plain_r = plain_l ^ (sl_2(plain_r, a) & sl_2(plain_r, b)) ^ sl_2(plain_r, c) ^ key

    return  plain_r





def dec(ciper_l, ciper_r, roundkey, round,a,b,c):
    for i in range(round):
        tmp=ciper_r
        ciper_r = onedec(ciper_l, ciper_r, roundkey[:,round-1-i],a,b,c)
        ciper_l=tmp
    return ciper_l, ciper_r


def array_to_int(a):
    b = 0
    for i in range(len(a)):
        b += a[i] << i
    return b

def int_to_array(a,block_n):
    arr=np.zeros(block_n, dtype=bool)
    for i in range(block_n):
        arr[i]=(a>>i)%2
    return arr


def veirfy_1(round, diff_in_l, diff_in_r, diff_out_l, diff_out_r,diff_mainkey_1,diff_mainkey_2,block_n,num_of_key,a,b,c,z):

    sum = 0

    for i in range(date_num//date_block):

        mainkey_1 = np.random.randint(0, 2, size=(date_block,num_of_key, block_n))
        mainkey_2 = sl_3(mainkey_1,offset)[:] ^ diff_mainkey_1

        roundkey_1 = keyschedule(mainkey_1, round,block_n,num_of_key,a,b,c,z)
        roundkey_2 = keyschedule(mainkey_2, round, block_n, num_of_key, a, b, c, z)


        ni_mainkey_3 = roundkey_1[:,round-num_of_key:] ^ diff_mainkey_2
        ni_mainkey_4 = roundkey_2[:,round - num_of_key:] ^ sl_2(diff_mainkey_2,offset)

        roundkey_3 = ni_keyschedule(ni_mainkey_3,round)
        roundkey_4 = ni_keyschedule(ni_mainkey_4,round)

        plain_l_1 = np.random.randint(0, 2, size=(date_block,block_n))
        plain_r_1 = np.random.randint(0, 2, size=(date_block,block_n))

        plain_l_2 = sl_2(plain_l_1,offset) ^ diff_in_l
        plain_r_2 = sl_2(plain_r_1,offset) ^ diff_in_r

        ciper_l_1, ciper_r_1 = enc(plain_l_1, plain_r_1, roundkey_1, round,a,b,c)
        ciper_l_2, ciper_r_2 = enc(plain_l_2, plain_r_2, roundkey_2, round,a,b,c)

        ciper_l_3 = ciper_l_1 ^ diff_out_l
        ciper_r_3 = ciper_r_1 ^ diff_out_r
        ciper_l_4 = ciper_l_2 ^ sl(diff_out_l,offset)
        ciper_r_4 = ciper_r_2 ^ sl(diff_out_r,offset)

        plain_l_3, plain_r_3 = dec(ciper_l_3,ciper_r_3,roundkey_3, round,a,b,c)
        plain_l_4, plain_r_4 = dec(ciper_l_4, ciper_r_4, roundkey_4, round,a,b,c)




        sum+=np.sum(~((sl_2(plain_l_3,offset) ^ plain_l_4^diff_in_l).any(axis=1)|(sl_2(plain_r_3,offset) ^ plain_r_4 ^ diff_in_r).any(axis=1)))
    print(sum)
    print(date_num)
    if sum!=0:
        print('weight of probility',math.log2(date_num / sum))
    return None




if __name__ == '__main__':





    diff_in_l_arr = int_to_array(diff_in_l_int,block_n)
    diff_in_r_arr = int_to_array(diff_in_r_int,block_n)
    diff_out_l_arr = int_to_array(diff_out_l_int,block_n)
    diff_out_r_arr = int_to_array(diff_out_r_int,block_n)
    diff_mainkey_1_arr = np.array([int_to_array(diff_mainkey_1[i],block_n)for i in range(num_of_key)])
    diff_mainkey_2_arr = np.array([int_to_array(diff_mainkey_2[i],block_n) for i in range(num_of_key)])

    veirfy_1(round, diff_in_l_arr, diff_in_r_arr, diff_out_l_arr, diff_out_r_arr,diff_mainkey_1_arr,diff_mainkey_2_arr,block_n,num_of_key,a,b,c,z[z_choose])


