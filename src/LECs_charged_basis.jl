




function LECs_cBorn_s(i::Int64, j::Int64)
    l = fill(0e0, 6)
    D, F = 0.8, 0.46
    if (i, j) == (1, 1)
        l[3] = (-D -3e0F)^2
        l[5] = (sqrt(3.)*(D-F) )^2
        return l
    elseif (i, j) == (1, 2) || (i, j) == (2, 1)
        l[3] = (-D - 3e0F) *(2e0D)
        l[5] = sqrt(3.)*(D-F) * (-2e0sqrt(3.)*F)
        return l
    elseif (i, j) == (1, 3) || (i, j) == (3, 1)
        l[3] = (-D - 3e0F) * (2e0D)
        l[5] = sqrt(3.) * (D - F) * 0
        return l
    elseif (i, j) == (1, 4) || (i, j) == (4, 1)
        l[3] = (-D - 3e0F) * (-D - 3e0F)
        l[5] = sqrt(3.) * (D - F) * (-sqrt(3.) * (D - F) )
        return l
    elseif (i, j) == (1, 5) || (i, j) == (5, 1)
        l[3] = (-D -3e0F) * (2e0D)
        l[5] = sqrt(3.)*(D-F) * (2e0sqrt(3.)*F)
        return l
    elseif (i, j) == (1, 6) || (i, j) == (6, 1)
        l[3] = (-D -3e0F) * 0
        l[5] = sqrt(3.) * (D - F) * (2e0D)
        return l
    elseif (i, j) == (2, 2)
        l[3] = 2e0D * 2e0D
        l[5] = -2e0sqrt(3.)F * (-2e0)sqrt(3.)F
        return l
    elseif (i, j) == (2, 3) || (i, j) == (3, 2)
        l[3] = 2e0D * 2e0D
        l[5] = -2e0sqrt(3.)F * 0
        return l
    elseif (i, j) == (2, 4) || (i, j) == (4, 2)
        l[3] = 2e0D * (-D -3e0F)
        l[5] = -2e0sqrt(3.)F * (-sqrt(3.)*(D-F) )
        return l
    elseif (i, j) == (2, 5) || (i, j) == (5, 2)
        l[3] = 2e0D * (2e0D)
        l[5] =-2e0sqrt(3.)F * 2e0sqrt(3)*F
        return l
    elseif (i, j) == (2, 6) || (i, j) == (6, 2)
        l[3] = 2e0D * 0
        l[5] = -2e0sqrt(3.)F * 2e0D
        return l
    elseif (i, j) == (3, 3)
        l[3] = 2e0D * 2e0*D
        l[5] = 0
        return l
    elseif (i, j) == (3, 4) || (i, j) == (4, 3)
        l[3] = 2e0D * (-D - 3e0F)
        l[5] = 0 * (-sqrt(3) * (D - F) )
        return l
    elseif (i, j) == (3, 5) || (i, j) == (5, 3)
        l[3] = 2e0D * 2e0D
        l[5] = 0 * (2sqrt(3)F)
        return l
    elseif (i, j) == (3, 6) || (i, j) == (6, 3)
        l[3] = 2e0D * 0
        l[5] = 0 * 2e0D
        return l
    elseif (i, j) == (4, 4) 
        l[3] = (-D - 3e0F)^2
        l[5] = (-sqrt(3.)* (D- F))^2
        return l
    elseif (i, j) == (4, 5) || (i, j) == (5, 4)
        l[3] = (-D - 3e0F) * 2e0D
        l[5] = (-sqrt(3.)* (D- F)) * 2e0sqrt(3.)F
        return l
    elseif (i, j) == (4, 6) || (i, j) == (6, 4)
        l[3] = (-D -3e0F) * 0
        l[5] = (-sqrt(3)*(D -F) ) * 2e0D
        return l
    elseif (i, j) == (5, 5)
        l[3] = (2e0D)^2
        l[5] = (2e0sqrt(3.)*F)^2
        return l
    elseif (i, j) == (5, 6) || (i, j) == (6, 5)
        l[3] = (2e0D) * 0
        l[5] = (2e0sqrt(3)F) * 2e0D
        return l
    elseif (i, j) == (6, 6)
        l[3] = 0
        l[5] = (2e0D)^2
        return l
    end
end

function LECs_cBorn_u(i::Int64, j::Int64)
    l = fill(0e0, 6)
    D, F = 0.8, 0.46
    if (i, j) == (1, 1)
        return l

    elseif (i, j) == (1, 2) || (i, j) == (2, 1)
        return l

    elseif (i, j) == (1, 3) || (i, j) == (3, 1)
        l[1] = sqrt(3.) * (D + F) * sqrt(3.) * (D - F)
        return l
    elseif (i, j) == (1, 4) || (i, j) == (4, 1)
        return l
    elseif (i, j) == (1, 5) || (i, j) == (5, 1)
        l[2] = sqrt(6.) * (D + F) * (sqrt(6.) * (D - F) )
        return l
    elseif (i, j) == (1, 6) || (i, j) == (6, 1)
        l[1] = sqrt(3.) * (D + F) * (-D -3e0F)
        return l
    elseif (i, j) == (2, 2)
        return l
    elseif (i, j) == (2, 3) || (i, j) == (3, 2)
        return l
    elseif (i, j) == (2, 4) || (i, j) == (4, 2)
        l[1] = sqrt(6.)*(D + F) * sqrt(6.) * (D - F)
        return l
    elseif (i, j) == (2, 5) || (i, j) == (5, 2)
        l[3] = (2e0D)^2
        return l
    elseif (i, j) == (2, 6) || (i, j) == (6, 2)
        return l
    elseif (i, j ) == (3, 3) 
        l[3] = 2e0D * 2e0D
        return l
    elseif (i, j) == (3, 4) || (i, j) == (4, 3)
        l[2] = -sqrt(3.) *(D - F) * (-sqrt(3.) * (D + F))
        return l
    elseif (i, j) == (3, 5) || (i, j) == (5, 3)
        return l
    elseif (i, j) == (3, 6) || (i, j) == (6, 3)
        return l

    elseif (i, j) == (4, 4)
        return l
    elseif (i,  j) == (4, 5) || (i, j) == (5, 4)
        return l
    elseif (i, j) == (4, 6) || (i, j) == (6, 4)
        l[2] = (-sqrt(3.)*(D+F)) * (-D-3e0F)
        return l
    elseif (i, j) == (5, 5)
        return l

    elseif (i, j) == (5, 6) || (i, j) == (6, 5)
        return l
    elseif (i, j) == (6, 6)
        l[5] = (2e0D)^2
        return l
    end
    
end