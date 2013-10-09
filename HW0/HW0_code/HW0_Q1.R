r = 1:100                         #Defining a vectore whose elements are 1 to 100
r[seq.int(3,99,3)]="Fizz"         #Assigning Fizz to multiples of 3 
r[seq.int(5,100,5)]="Buzz"        #Assigning Buzz to multiples of 5 
r[seq.int(15,90,15)]="FizzBuzz"   #Assigning FizzBuzz to multiples of 15
print(r)