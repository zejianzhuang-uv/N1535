

function kcm(w, m1, m2)
    if w > m1 + m2
        return sqrt(kallen(w^2, m1^2, m2^2) ) / (2w)
    else
        return 0e0
    end
end


function phase_space(w, par, ch)
    s = w * w
    mch = par[:mch]
    M1 = mch[1][2]
    kcm1 = kcm(w, mch[1]...)
    if ch == :ch11
        psp = 1 / (4π) * M1 * mch[1][2] / s * (kcm(w, mch[1]...) / kcm1)
    elseif ch == :ch12
        psp = 1 / (4π) * M1 * mch[2][2] / s * kcm(w, mch[2]...) / kcm1
        
    elseif ch == :ch13
        psp = 1 / (4π) * M1 * mch[3][2] / s * kcm(w, mch[3]...) / kcm1
        
    elseif ch == :ch14

        psp = 1 / (4π) * M1 * mch[4][2] / s * kcm(w, mch[4]...) / kcm1
        
    elseif ch == :ch15
        psp = 1 / (4π) * M1 * mch[5][2] / s * kcm(w, mch[5]...) / kcm1
        
    elseif ch == :ch16
        psp = 1 / (4π) * M1 * mch[6][2] / s * kcm(w, mch[6]...) / kcm1
        
    end
    # In a physical process, the momentum should be greater than 0
    if w > sum(mch[1])
        return psp
    else
        return 0e0
    end
end