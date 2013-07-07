! Melting pressure [bars] vs T [mK]
! Arg: T = 0 .. 250 [mK]
! Note: Wrong for T<0.2mK?? - check range!
! Ref: Greywall. Phys. Rev.B v.33 #11 p.7520
! Origin: Mukharskii, Dmitriev

      function He3_Pmelt(T)
        implicit none
        include '../he3.fh'
        real*8 T
        if (T.gt.0.and.T.le.250) then
          He3_Pmelt = He3_Pa
     .     - 0.19652970D-1*T**(-3)
     .     - 0.61880268D-1*T**(-2)
     .     - 0.78803055D-1*T**(-1)
     .     + 0.13050600D0
     .     - 0.43519381D-1*T
     .     + 0.13752791D-3*T**2
     .     - 0.17180436D-6*T**3
     .     - 0.22093906D-9*T**4
     .     + 0.85450245D-12*T**5
        else
          He3_Pmelt = -1D0
        endif
      end
