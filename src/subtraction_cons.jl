


function subtraction_cons(mchi, mu, qmax)
    m1, m2 = mchi
    a = -2/(m1+m2) * ( m1*log(1 + sqrt(1 + m1^2/qmax^2)) + m2*log(1 + sqrt(1 + m2^2/qmax^2) ) ) + 2log(mu / qmax)
    return a
end

function subtraction_cons_tuple(params; mu=0, qmax=0, n=2)
    if mu == 0 || qmax == 0
        return 0
    else
        return Tuple(subtraction_cons(mchi, mu, qmax) for mchi in params[:mch] )
    end
end



