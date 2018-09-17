      function uxx=dss042(xl,xu,n,u,ux,nl,nu)
%
%   SUBROUTINE DSS042 COMPUTES A SECOND-ORDER APPROXIMATION OF A
%   SECOND-ORDER DERIVATIVE, WITH OR WITHOUT THE NORMAL DERIVATIVE
%   AT THE BOUNDARY.
%
%   ARGUMENT LIST
%
%      XL      LEFT VALUE OF THE SPATIAL INDEPENDENT VARIABLE (INPUT)
%
%      XU      RIGHT VALUE OF THE SPATIAL INDEPENDENT VARIABLE (INPUT)
%
%      N       NUMBER OF SPATIAL GRID POINTS, INCLUDING THE END
%              POINTS (INPUT)
%
%      U       ONE-DIMENSIONAL ARRAY OF THE DEPENDENT VARIABLE TO BE
%              DIFFERENTIATED (INPUT)
%
%      UX      ONE-DIMENSIONAL ARRAY OF THE FIRST DERIVATIVE OF U.
%              THE END VALUES OF UX, UX(1) AND UX(N), ARE USED IN
%              NEUMANN BOUNDARY CONDITIONS AT X = XL AND X = XU,
%              DEPENDING ON THE ARGUMENTS NL AND NU (SEE THE DE-
%              SCRIPTION OF NL AND NU BELOW)
%
%      UXX     ONE-DIMENSIONAL ARRAY OF THE SECOND DERIVATIVE OF U
%              (OUTPUT)
%
%      NL      INTEGER INDEX FOR THE TYPE OF BOUNDARY CONDITION AT
%              X = XL (INPUT).  THE ALLOWABLE VALUES ARE
%
%                 1 - DIRICHLET BOUNDARY CONDITION AT X = XL
%                     (UX(1) IS NOT USED)
%
%                 2 - NEUMANN BOUNDARY CONDITION AT X = XL
%                     (UX(1) IS USED)
%
%      NU      INTEGER INDEX FOR THE TYPE OF BOUNDARY CONDITION AT
%              X = XU (INPUT).  THE ALLOWABLE VALUES ARE
%
%                 1 - DIRICHLET BOUNDARY CONDITION AT X = XU
%                     (UX(N) IS NOT USED)
%
%                 2 - NEUMANN BOUNDARY CONDITION AT X = XU
%                     (UX(N) IS USED)
%
%
%   THE FOLLOWING DERIVATION IS FOR A SET OF SECOND-ORDER, FOUR-POINT
%   APPROXIMATIONS FOR A SECOND DERIVATIVE THAT CAN BE USED AT THE
%   BOUNDARIES OF A SPATIAL DOMAIN.  THESE APPROXIMATIONS HAVE THE
%   FEATURES
%
%      (1)  ONLY INTERIOR AND BOUNDARY POINTS ARE USED (I.E., NO
%           FICTITIOUS POINTS ARE USED)
%
%      (2)  THE NORMAL DERIVATIVE AT THE BOUNDARY IS INCLUDED AS PART
%           OF THE APPROXIMATION FOR THE SECOND DERIVATIVE
%
%      (3)  APPROXIMATIONS FOR THE BOUNDARY CONDITIONS ARE NOT USED.
%
%   THE DERIVATION IS BY PROFESSOR GILBERT A. STENGLE, DEPARTMENT OF
%   MATHEMATICS, LEHIGH UNIVERSITY, BETHLEHEM, PA 18015, AND WAS DONE
%   ON DECEMBER 7, 1985.
%
%   FOR AN APPROXIMATION AT THE LEFT BOUNDARY, INVOLVING THE POINTS
%   I, I+1, I+2 AND I+3, CONSIDER THE FOLLOWING TAYLOR SERIES EXPAN-
%   SIONS
%
%                   UX(I)( DX)   UXX(I)( DX)**2   UXXX(I)( DX)**3
%   U(I+1) = U(I) + ---------- + -------------- + --------------- +...
%                        1             2                 6
%
%
%                   UX(I)(2DX)   UXX(I)(2DX)**2   UXXX(I)(2DX)**3
%   U(I+2) = U(I) + ---------- + -------------- + --------------- +...
%                        1             2                 6
%
%   IF WE NOW FORM THE FOLLOWING LINEAR COMBINATION, INVOLVING CON-
%   STANTS A, B, C AND D TO BE DETERMINED, AND USE THE PRECEDING TWO
%   TAYLOR SERIES,
%
%      A*U(I) + B*UX(I) + C*U(I+1) + D*U(I+2)
%
%   WE HAVE
%
%      A*U(I) + B*UX(I) + C*U(I+1) + D*U(I+2) =
%
%      (A + B + C + D)*U(I) +
%
%      (B + DX*C + 2*DX*D)*UX(I) +
%
%      (C*(DX**2)/2 + D*((2*DX)**2)/2)*UXX(I) +
%
%      (C*(DX**3)/6 + D*((2*DX)**3)/6)*UXXX(I) + O(DX**4)
%
%   THE THIRD DERIVATIVE, UXXX(I), CAN BE DROPPED BY TAKING
%
%      C = -8*D
%
%   THE SECOND DERIVATIVE, UXX(I), CAN BE RETAINED BY TAKING
%
%      (DX**2)(C/2 + 2*D) = 1
%
%   WHICH, WHEN COMBINED WITH THE PRECEDING RESULT GIVES
%
%      D = -1/(2*(DX**2))
%
%      C = 4/(DX**2)
%
%   THE FIRST DERIVATIVE, UX(I), CAN BE DROPPED BY TAKING
%
%      B + DX*C + 2*DX*D = 0
%
%   OR
%
%      B = -DX*C - 2*DX*D = -4/DX - 2*DX*(-1/(2*(DX**2))) = -3/DX
%
%   FINALLY, U(I), CAN BE DROPPED BY TAKING
%
%      A = - C - D = 8*D - D = -7*D = -7/(2*(DX**2))
%
%   IF WE NOW SOLVE FOR THE DERIVATIVE OF INTEREST, UXX(I),
%
%      UXX(I) = -7/(2(DX**2))*U(I) - 3/DX*UX(I)
%
%               + 8/(DX**2)*U(I+1) - 1/(2*(DX**2))U(I+2) + O(DX**2)
%
%        = (1/(2*(DX**2)))*(-U(I+2) + 8*U(I+1) - 7*U(I) - 6*DX*UX(I))
%
%          + O(DX**2)
%
%   WHICH IS THE FOUR-POINT, SECOND-ORDER APPROXIMATION FOR THE SECOND
%   DERIVATIVE, UXX(I), INCLUDING THE FIRST DERIVATIVE, UX(I).
%
%   FOUR CHECKS OF THIS APPROXIMATION CAN EASILY BE MADE FOR U(I) =
%   1, U(I) = X, U(I) = X**2 AND U(I) = X**3
%
%      UXX(I) = (1/(2*(DX**2)))*(-1 + 8*1 - 7*1 - 6*DX*0) = 0
%
%      UXX(I) = (1/(2*(DX**2)))*(-(X + 2*DX) + 8*(X + DX)
%
%               -7*X - 6*DX*1) = 0
%
%      UXX(I) = (1/(2*(DX**2)))*(-(X + 2*DX)**2 + 8*(X + DX)**2
%
%             - 7*(X**2) - 6*DX*(2*X))
%
%              = (-  X**2 -  4*X*DX - 4*DX**2
%
%                + 8*X**2 + 16*X*DX + 8*DX**2
%
%                - 7*X**2 - 12*X*DX)/(2*(DX**2)) = 2
%
%      UXX(I) = (1/(2*(DX**2)))*(-(X + 2*DX)**3 + 8*(X + DX)**3
%
%             - 7*(X**3) - 6*DX*(3*X**2))
%
%             = (1/(2*(DX**2)))*(- X**3 - 6*DX*X**2 - 12*X*DX**2
%
%             - 8*DX**3 + 8*X**3 + 24*DX*X**2 + 24*X*DX**2 + 8*DX**3
%
%             - 7*X**3 - 18*DX*X**2)
%
%             = (1/(2*(DX**2)))*(12*X*DX**2) = 6*X
%
%   THE PRECEDING APPROXIMATION FOR UXX(I) CAN BE APPLIED AT THE
%   LEFT BOUNDARY VALUE OF X BY TAKING I = 1.  AN APPROXIMATION AT
%   THE RIGHT BOUNDARY IS OBTAINED BY TAKING DX = -DX AND REVERSING
%   THE SUBSCRIPTS IN THE PRECEDING APPROXIMATION, WITH I = N
%
%      UXX(I)
%
%        = (1/(2*(DX**2)))*(-U(I-2) + 8*U(I-1) - 7*U(I) + 6*DX*UX(I))
%
%          + O(DX**2)
%
%   TO OBTAIN APPROXIMATIONS OF THE SECOND DERVIAVTIVE WHICH DO NOT
%   INVOLVE THE FIRST DERIVATIVE, WE TAKE AS THE LINEAR COMBINATION
%
%      A*U(I) + B*U(I+1) + C*U(I+2) + D*U(I+3) =
%
%   WE HAVE
%
%      A*U(I) + B*U(I+1) + C*U(I+2) + D*U(I+3) =
%
%      (A + B + C + D)*U(I)+
%
%      (DX*B + 2*DX*C + 4*DX*D)*UX(I)+
%
%      (B*(DX**2)/2 + C*((2*DX)**2)/2 + D*((3*DX)**2)/2)*UXX(I) +
%
%      (B*(DX**3)/6 + C*((2*DX)**3)/6 + D*((3*DX)**3)/6)*UXX(I) +
%
%      O(DX**4)
%
%   THE THIRD DERIVATIVE, UXXX(I), CAN BE DROPPED BY TAKING
%
%      B + 8*C + 27*D = 0
%
%   THE SECOND DERIVATIVE, UXX(I), CAN BE RETAINED BY TAKING
%
%      (DX**2)*(B/2 + 2*C + (9/2)*D) = 1
%
%   THE FIRST DERIVATIVE CAN BE DROPPED BY TAKING
%
%      B + 2*C + 3*D = 0
%
%   SOLUTION OF THE PRECEDING EQUATIONS FOR C AND D BY ELIMINATION OF
%   B GIVES
%
%      6*C + 24*D = 0
%
%      4*C + 18*D = -2/(DX**2)
%
%   THEN, ELIMINATING C, GIVES
%
%      (18 - 16)*D = -2/(DX**2)
%
%   OR
%
%      D = -1/(DX**2)
%
%      C = (24/6)/(DX**2) = 4/(DX**2)
%
%      B = -8/(DX**2) + 3/(DX**2) = -5/(DX**2)
%
%   U(I) CAN BE DROPPED BY TAKING
%
%      A + B + C + D = 0
%
%   OR
%
%      A = (5 - 4 + 1)/(DX**2) = 2/(DX**2)
%
%   IF WE NOW SOLVE FOR THE DERIVATIVE OF INTEREST, UXX(I),
%
%      UXX(I) = (1/DX**2)*(2*U(I) - 5*U(I+1) + 4*U(I+2) - 1*U(I+3))
%
%             + O(DX**2)
%
%   WHICH IS THE FOUR-POINT, SECOND-ORDER APPROXIMATION FOR THE SECOND
%   DERIVATIVE, UXX(I), WITHOUT THE FIRST DERIVATIVE, UX(I).
%
%   FOUR CHECKS OF THIS APPROXIMATION CAN EASILY BE MADE FOR U(I) =
%   1, U(I) = X, U(I) = X**2 AND U(I) = X**3
%
%      UXX(I) = (1/DX**2)*(2 - 5 + 4 - 1) = 0
%
%      UXX(I) = (1/DX**2)*(2*X - 5*(X + DX) + 4*(X + 2*DX)
%
%               - 1*(X + 3*DX)) = 0
%
%      UXX(I) = (1/DX**2)*(2*X**2 - 5*(X + DX)**2 + 4*(X + 2*DX)**2
%
%               - 1*(X + 3*DX)**2) = 2
%
%      UXX(I) = (1/DX**2)*(2*X**3 - 5*(X + DX)**3 + 4*(X + 2*DX)**3
%
%             - 1*(X + 3*DX)**3)
%
%              = (1/DX**2)*(2*X**3 - 5*X**3 - 15*X*DX**2
%
%              - 15*DX*X**2 - 5*DX**3 + 4*X**3 + 24*DX*X**2
%
%              + 48*X*DX**2 + 32*DX**3 - X**3 - 9*DX*X**2
%
%              - 27*X*DX**2 - 27DX**3)
%
%              = (1/DX**2)*(6*X*DX**2) = 6*X
%
%   THE PRECEDING APPROXIMATION FOR UXX(I) CAN BE APPLIED AT THE
%   LEFT BOUNDARY VALUE OF X BY TAKING I = 1.  AN APPROXIMATION AT
%   THE RIGHT BOUNDARY IS OBTAINED BY TAKING DX = -DX AND REVERSING
%   THE SUBSCRIPTS IN THE PRECEDING APPROXIMATION, WITH I = N
%
%      UXX(I) = (1/DX^2)*(2*U(I) - 5*U(I-1) + 4*U(I-2) - 1*U(I-3))
%
%             + O(DX^2)
%
%   GRID SPACING
      dx=(xu-xl)/(n-1);
%
%   calculate uxx at the left boundary, without ux
      if(nl==1),
         uxx(1)=(2*u(1)-5*u(2)+4*u(3)-u(4))/(dx^2);
%
%   calculate uxx at the left boundary, including ux
      else
         uxx(1)=(-7*u(1)+8*u(2)-u(3))/(2*dx^2)-6*ux(1)/(2*dx);
      end
%
%   calculate uxx at the right boundary, without ux
      if(nu==1),
         uxx(n)=(2*u(n)-5*u(n-1)+4*u(n-2)-u(n-3))/(dx^2);
%
%   calculate uxx at the right boundary, including ux
      else
         uxx(n)=(-7*u(n)+8*u(n-1)-u(n-2))/(2*dx^2)+6*ux(n)/(2*dx);
      end
%
%   calculate uxx at the interior grid points
      uxx(2:n-1)=(u(3:n)-2*u(2:n-1)+u(1:n-2))/dx^2;
     







