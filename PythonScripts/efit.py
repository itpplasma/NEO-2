class efit_type:
  characters_per_number=16
  numbers_per_line=5

  def __init__(self, filename):
    import numpy as np

    # Read in the file.
    with open(filename) as f:
        lines = f.readlines()

    # Split the first line, which contains some usefull (e.g. nw, nh) and
    # some not so usefull information (e.g. unknownvaluestr and signtodiscard).
    self.unknownvaluestr, self.datestr, self.signtodiscard, self.shootnr, self.timestr, self.idum,self.nw,self.nh = lines[0].split()

    # Variables are so far all strings, but change to integer for some of them.
    self.shootnr = int(self.shootnr)
    self.idum = int(self.idum)# Some
    self.nw = int(self.nw)# Number of horizontal R grid points
    self.nh = int(self.nh)# Number of vertical z grid points

    # First few (four) lines contain specific values.
    self.rdim, self.zdim, self.rcentr, self.rleft, self.zmid = np.array(lines[1].split(),dtype=float)
    self.rmaxis, self.zmaxis, self.simag1, self.sibry1, self.bcentr = np.array(lines[2].split(),dtype=float)
    self.current, self.simag2, self.dummy1, self.rmaxis, self.dummy2 = np.array(lines[3].split(),dtype=float)
    self.zmaxis, self.dummy3, self.sibry2, self.dummy4, self.dummy5 = np.array(lines[4].split(),dtype=float)
    #~ self.rmaxis.unit = "m"
    #~ self.zmaxis.unit = "m"
    #~ self.simag.unit = "Weber/rad"
    #~ self.sibry.unit = "Weber/rad"

    # Create array of psi points.
    self.psi = np.linspace(0,self.sibry1-self.simag1,self.nw)

    # Now read in six array, the first and the last one of size nw, the fourth of size nw*nh.
    offset = 5
    self.fpol = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))

    offset = offset+int(self.nw/self.numbers_per_line)
    self.pres = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))
    #~ self.pres.unit = "N/m^2"

    offset = offset+int(self.nw/self.numbers_per_line)
    self.ffprim = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))
    #~ self.ffprim.unit = "(mT)^2/(Weber/rad)"

    offset = offset+int(self.nw/self.numbers_per_line)
    self.pprime = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))
    #~ self.pprime.unit = "(N/m^2)/(Weber/rad)"

    offset = offset+int(self.nw/self.numbers_per_line)
    self.psi2d = self.local_split(''.join(lines[offset:offset+int(self.nw*self.nh/self.numbers_per_line)])).reshape([self.nh,self.nw])
    #~ self.psi2d.unit = "Weber/rad"

    offset = offset+int(self.nw*self.nh/self.numbers_per_line)
    self.q = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))

    # Next are two sizes(?) in the file.
    offset = offset+int(self.nw/self.numbers_per_line)
    self.number_boundary_points, self.number_limiter_points = lines[offset].split()
    self.number_boundary_points = int(self.number_boundary_points)
    self.number_limiter_points = int(self.number_limiter_points)

    # It follows an array of (r,z) values of the boundary points.
    offset = offset+1
    self.unknown_array1 = self.local_split(''.join(lines[offset:offset+int(np.ceil(2*self.number_boundary_points/self.numbers_per_line))])).reshape([self.number_boundary_points,2])

    # The following code does not work, as the format of the file is not
    # as expected.
    # For the default file and another one there are 4*5 additional
    # values (0) at the end of the array.
    #~ offset = offset+int(2*self.number_boundary_points/self.numbers_per_line)
    #~ self.unknown_array1 = self.local_split(''.join(lines[offset:offset+int(np.ceil(2*self.number_limiter_points/self.numbers_per_line))])).reshape([self.number_limiter_points,2])

    # Rest of the the file has so far unknown content. There seems to
    # follow a line with three values, another array and
    # a last value?
    # If it is known what this is, it is read, until then we just ignore
    # it.

  def local_split(self, a):
    import numpy as np

    aflat = a.replace('\n','')
    alist = [aflat[self.characters_per_number*k:self.characters_per_number*(k+1)] for k in range(int(len(aflat)/self.characters_per_number))]
    return np.array(alist, dtype=float)

########################################################################
if __name__ == "__main__":
  import matplotlib.pyplot as plt

  # If this is run as main, just create a class-object with default
  # name(s) for the file(s).

  filename = '/proj/plasma/Neo2/ASDEX-U/31114/g031114.05000'
  #~ filename = '/proj/plasma/Neo2/ASDEX-U/32566/g032566.05000'
  efitvar = efit_type(filename)

  plt.figure()
  plt.plot(efitvar.psi,efitvar.fpol)
  plt.xlabel('psi_pol')
  plt.ylabel('F_pol')
  plt.savefig('F_pol')

  plt.figure()
  plt.plot(efitvar.psi,efitvar.pres)
  plt.xlabel('psi_pol')
  plt.ylabel('pressure')
  plt.savefig('pressure')

  plt.figure()
  plt.plot(efitvar.psi,efitvar.ffprim)
  plt.xlabel('psi_pol')
  plt.ylabel('F*dF/dpsi_pol')
  plt.savefig('FdF_dpsi_pol')

  plt.figure()
  plt.plot(efitvar.psi,efitvar.pprime)
  plt.xlabel('psi_pol')
  plt.ylabel('dp/dpsi_pol')
  plt.savefig('dp_dpsi_pol')

  plt.figure()
  plt.contour(efitvar.psi2d,30)
  plt.colorbar()
  plt.axis('equal')
  plt.savefig('equal')

  plt.figure()
  plt.plot(efitvar.psi,efitvar.q)
  plt.ylim(0,5)
  plt.xlabel('psi_pol')
  plt.ylabel('safety factor q')
  plt.savefig('safety_factor')
