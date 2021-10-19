for A in range(1,10000):
    count_2= 0
    for x in range (1,10000):
        if x%A!=0 or x%16==0 or x%23==0:
            count_2+=1
    if count_2 == 9999:
        print(A)
        
print("end")