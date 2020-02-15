


       do J = 1, NUPAR
          CUPAR=UPAR(J)

	
	    g_33j = 0.0

          do K = 2, NUPER - 1
             CUPER=UPER(K)
             g_33j = g_33j + UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))
          end do

	    k = 1
          g_33j = g_33j + 0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))

	    k = nuper
          g_33j = g_33j + 0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))
	
	    du = (uper(nuper) - uper(1))/(nuper - 1)

	    ga_33(j) = g_33j * du

       end do

!      the uperp integral is now done; now do parallel


       call EQSIMPSON1D_2(NUPAR, UPAR, ga_33, g_33)
	



