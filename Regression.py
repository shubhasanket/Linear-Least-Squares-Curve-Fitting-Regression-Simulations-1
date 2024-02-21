"""Linear Least Squares Curve Fitting (Regression) Simulations
"""
import numpy as n
import xlrd
import matplotlib.pyplot as plt

class LeastSquaresCurveFitting:
     def __init__ (self, data, y_set = 0, dim = 3):
          """data: a 2-d list:
             data[0] contains the x values;
             data[1] is the list which contains the y values.
             y_set: the value of the set y in use, (0-3)
             dim: (dim - 1) is the order of the polynomial
          """

          self.y_set = y_set
          self.data = n.array(data)
          self.dim = dim
          self.A = self.create_A()# This is the Vandermonde matrix 
          self.x_hat = self.compute_x_hat() # Least squares solutions x hat
          
          # e is the error vector calculated with the formula b-Ax
          # TLSE stands for Total Least Squares Error
          self.e_bAx, self.TLSE = self.compute_e_and_TLSE() 

          # Projection matrix P
          self.P = self.compute_P()

          # matrix I-P
          self.I_minus_P = self.compute_I_minus_P()
          
          #vector e calculated with the equation (I-P)b 
          self.e_IPb = n.array(n.dot(self.I_minus_P, self.data[1]))

          self.equation = ''

          #This generates the equation of the line
          self.plot(re=True)

          #These 2 lists will be used while plotting the graph
          self.x_axis_vals = []
          self.y_axis_vals = []
          
     def create_A(self):
          """This method creates the Vandermonde matrix
          """
          l = []
          for i in range (len(self.data[0])):
               l.append([])
               for j in range (self.dim):
                    l[i].append(self.data[0][i]**j)
          return n.array(l)


     def compute_x_hat(self):
          """This method calculates x_hat
             x_hat = ((A_T * A)^-1)*A_T*b
             A_T is the transpose of A
          """
          a = n.dot(self.A.T, self.A) # A_T * A
          b = n.linalg.inv(a) # (A_T * A)^-1
          d = n.dot(b, n.dot(self.A.T, self.data[1]))
          return d

          
     def compute_P(self):
          """This method calculates the projection matrix P
             P = A*((A_T * A)^-1)*A_T
          """
          a = n.linalg.inv(n.dot(self.A.T, self.A)) #(A_T * A)^-1
          return n.dot(n.dot(self.A, a), self.A.T) #A((A_T * A)^-1)A_T


     def compute_e_and_TLSE(self, l=[]):
          """This method calculates the error vector e and TLSE
          """
          if l == []:
               b = self.data[1]
          else:
               b = l
          e = n.subtract(b, n.dot(self.A,self.x_hat))#e = b-Ax
          return e, n.dot(e,e)


     def compute_I_minus_P(self):
          """This method calculates I-P
             This matrix can also be used to calculate e
          """
          i = n.identity(len(self.data[0]))
          return n.subtract(i,self.P)

     
     def plot(self, re=False, test_data = [], set_no = -1):
          """This method generates a graph for the best fit curve
             and returns the the equation of the line
          """
          #1st coefficient of the equation
          self.equation = str(format(self.x_hat[0],".2e"))
          #This loop generates the equation
          for i in range(1, len(self.x_hat)):
               
               if not self.x_hat[i]: # If a coefficient is 0
                    continue         # Skip
               
               else:
                    self.equation += " + " + str(format(self.x_hat[i],".2e"))\
                                     + "x^"+str(i)
                    
          if re:# If the function was invoked to calculate the equation of the line,
               return None # return None and exit
          """Kamal bhai, I know this return statement is inconsistent.
             As it is a conditional return, I would have to have another
             return somewhere below, but I did not know how to
             incorporate that. Could you please suggest something?
          """

          if test_data == []: # If no test_data has been provided, test_data = self.data
               test_data = self.data
          if set_no == -1: # This is to know which set of y's is being plotted
               set_no = self.y_set
          # plt.scatter plots the dots
          plt.scatter(n.array(test_data[0]), n.array(test_data[1]))

          # Kamal bhai, I do not know how the following command exactly works.
          # I saw similar statements in a few tutorials and examples and experimented
          # a bit to get to this.
          x = n.linspace(self.data[0][0], self.data[0][-1],
                         num = round((self.data[0][-1]-self.data[0][0])*25))

          y = 0
          for i in range (len(self.x_hat)):
               y += self.x_hat[i] * x ** i

          # Ploting the curve
          plt.plot(x, y, label = self.equation)

          #labels the axes
          plt.xlabel("x")
          plt.ylabel("y"+str(set_no))
          
          plt.title("Graph for set y"+str(set_no)+" (degree = "+str(self.dim)+")")
          plt.grid()
          print("Close the graph to continue")
          plt.show()
     
          
"""The following lines of code answer the questions in the PDF.
"""
DATA_SETS = 5
# Kamal bhai, you would have to change the path to run the program.
# It is the location of the excel file
LOC = r"C:\Users\Shubham\AppData\Local\Programs\Python\Python39\Workspace\Linear_Algebra\1\line_data.xls" # name of the excel file

def main():
     # Part 1
     data = read_excel(LOC) # data is a 2 dimensional list
     print("************** Part 1 **************\n")
     a = LeastSquaresCurveFitting(data = data[:2], y_set=0)
     print("The equation of the line is:"+a.equation)
     a.plot()
     print("Q: Is b = y0 in C(A)?\nA: Yes, it is.\nQ:How can you verify it?")
     print("b = A(x_hat) = \n", n.dot(a.A,a.x_hat), "\nand y0 =\n", a.data[1])
     print("Thus, b = y0.")
     print("TLSE = ", a.TLSE, "which can be rounded to", format(a.TLSE, ".1f"))
     print("The projection matrix P is:\n", a.P)
     print("P(transpose) =\n", a.P.T)
     print("P^2 =\n", n.dot(a.P,a.P))

     print("Therefore, P = P(transpose) = P^2")
     print("Matrix I-P =\n", a.I_minus_P)
     print("(I-P)(transpose) =\n", a.I_minus_P.T)
     print("(I-P)^2 =\n", n.dot(a.I_minus_P, a.I_minus_P))
     print("Therefore, I-P = (I-P)(transpose) = (I-P)^2")
     print("e = (I-P)b =\n", a.e_IPb)
     print("e = b-Ax = \n", a.e_bAx)
     print("The e(s) match. The minute discrepancies are caused by round-off errors.")
     print("As the individual components of e are 0, their average is 0 as well")
     print("e(transpose)*A =\n", n.dot(a.e_bAx, a.A))
     print("Thus, e lies in the nullspace of A")
     print("\n************** Part 2 **************\n")

     # Creating a list of instances of the class with differing parametres 
     l = []
     for i in range (2, DATA_SETS):
          l.append(LeastSquaresCurveFitting(data = [data[0],data[i]], y_set =(i-1), dim = 2))
          print("The following results are for the values y" +str(i-1))
          print("The equation of the line is:"+l[i-2].equation)
          l[i-2].plot()
          print("TLSE = ", l[i-2].TLSE, "which can be rounded to", format(l[i-2].TLSE, ".3f"))
          print("The projection matrix P is:\n", l[i-2].P)
          print("P(transpose) =\n", l[i-2].P.T)
          print("P^2 =\n", n.dot(l[i-2].P,l[i-2].P))
          print("Therefore, P = P(transpose) = P^2")
          print("Matrix I-P =\n", l[i-2].I_minus_P)
          print("(I-P)(transpose) =\n", l[i-2].I_minus_P.T)
          print("(I-P)^2 =\n", n.dot(l[i-2].I_minus_P, l[i-2].I_minus_P))
          print("Therefore, I-P = (I-P)(transpose) = (I-P)^2")
          print("e = (I-P)b =\n", l[i-2].e_IPb)
          print("e = b-Ax = \n", l[i-2].e_bAx)
          print("The e(s) match. The minute discrepancies are caused by round-off errors.")
          print("The average value of e is "+ format(sum(l[i-2].e_bAx)/len(l[i-2].e_bAx),".2f"))
          print("e(transpose)*A =\n", n.dot(a.e_bAx, a.A))
          print("Thus, e lies in the nullspace of A")
          print("Line parameters Vs true values")
          p = n.dot(l[i-2].A, l[i-2].x_hat)
          for j in range (len(p)):
               print(round(p[j],4), "|", round(l[i-2].data[1][j],4))
          print()

     print("\n************** Part 3 **************\n")
     print("Checking how well the model trained on the y2 data set fits the y3 data set(unseen/test data)")
     l[1].plot(test_data = [data[0], data[4]], set_no = 3)
     a,b = l[1].compute_e_and_TLSE(l = data[4])
     print("TLSE =", round(b,3), "and TLSE calculated in part 2 is", round(l[1].TLSE,3))
     print("Checking how well the model trained on the y3 data set fits the y2 data set(unseen/test data)")
     l[2].plot(test_data = [data[0], data[3]], set_no = 2)
     a,b = l[2].compute_e_and_TLSE(l = data[3])
     print("TLSE =", round(b,3), "and TLSE calculated in part 2 is", round(l[2].TLSE,3))
     print("Q: Are they very different?")
     print("A: Yes, they are very different")
     print("Q: Would you say that the models do about as well on the test data as on the training data?")
     print("A: No\n")
     print()

     print("************** Part 4 **************\n")
     
     d = [2,3,4,8,12,16]
     d = [i for i in range (2, 20)]
     lt = []
     #Creating a list of instances of the class with differing parametres
     for i in range (len(d)):
          lt.append([])
          for j in range(1, DATA_SETS):
               lt[i].append(LeastSquaresCurveFitting(data = [data[0],data[j]], y_set =(j-1), dim = d[i]))
          print("The equation of the line for data 0 and degree", d[i], "is:", lt[i][0].equation)
     print("To conclude, the equation of the line for data set y0 and for all degrees greater or equal to 2 is: ", end = '')
     print(lt[0][0].equation, "because x vs y0 is a perfect-fit straight line in 2 dimensions.\n")
     # Table of TLSE
     print("Table of TLSEs:")
     for i in range (DATA_SETS - 1):
          print("{0:^10}".format(" "*i+"y"+str(i)), end = '')
     print()
     
     for i in range (len(d)):
          print(d[i], end = '  ')
          print("{0:7} {1:^7} {2:^7} {3:^7}"
                .format(str(round(lt[i][0].TLSE,2)), str(round(lt[i][1].TLSE,2)), "   "+str(round(lt[i][2].TLSE,2)),
                        "   "+str(round(lt[i][3].TLSE,2))))
          print()

     #In the following section, one can choose any polynomial
     # degree that features in the list d and the data set to
     # be used for plotting.
     
     val = True    
     while val:
          x, y = eval(input("Enter the degree and the data set no.(e.g 12,3): "))
          if x == -1:#-1 to be entered to end the program
               val = False
          else:
               index = d.index(x)
               lt[index][y].plot()
          
     print("The program has ended")


def read_excel(loc):
     """Method to read the data from the excel sheet
     """
     # To open Workbook
     wb = xlrd.open_workbook(loc)
     sheet = wb.sheet_by_index(0)

     # Create lists by reading the data from the file
     x = [] 
     for i in range (DATA_SETS):
          x.append(sheet.col_values(i)[1:22])
          
     return x

main()










