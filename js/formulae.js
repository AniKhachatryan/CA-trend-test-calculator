function LogGamma(Z) {
    with (Math) {
        var S=1+76.18009173/Z-86.50532033/(Z+1)+24.01409822/(Z+2)-1.231739516/(Z+3)+.00120858003/(Z+4)-.00000536382/(Z+5);
        var LG= (Z-.5)*log(Z+4.5)-(Z+4.5)+log(S*2.50662827465);
    }
    return LG
}

function Gcf(X,A) {        // Good for X>A+1
    with (Math) {
        var A0=0;
        var B0=1;
        var A1=1;
        var B1=X;
        var AOLD=0;
        var N=0;
        while (abs((A1-AOLD)/A1)>.00001) {
            AOLD=A1;
            N=N+1;
            A0=A1+(N-A)*A0;
            B0=B1+(N-A)*B0;
            A1=X*A0+N*A1;
            B1=X*B0+N*B1;
            A0=A0/B1;
            B0=B0/B1;
            A1=A1/B1;
            B1=1;
        }
        var Prob=exp(A*log(X)-X-LogGamma(A))*A1;
    }
    return 1-Prob
}

function Gser(X,A) {        // Good for X<A+1.
    with (Math) {
        var T9=1/A;
        var G=T9;
        var I=1;
        while (T9>G*.00001) {
            T9=T9*X/(A+I);
            G=G+T9;
            I=I+1;
        }
        G=G*exp(A*log(X)-X-LogGamma(A));
    }
    return G
}

function Gammacdf(x,a) {
    var GI;
    if (x<=0) {
        GI=0
    } else if (x<a+1) {
        GI=Gser(x,a)
    } else {
        GI=Gcf(x,a)
    }
    return GI
}

function compute(x, df_val) {
    Z=eval(x)
    DF=eval(df_val)
    if (DF<=0) {
        alert("Degrees of freedom must be positive")
    } else {
        Chisqcdf=Gammacdf(Z/2,DF/2)
    }
    Chisqcdf = 1 - Math.round(Chisqcdf*100000)/100000;
    return Chisqcdf;
}

function Calculate() {
    var cases_AA = parseFloat(document.getElementById("case_AA").value);
    var cases_Aa = parseFloat(document.getElementById("case_Aa").value);
    var cases_aa = parseFloat(document.getElementById("case_aa").value);
    var controls_AA = parseFloat(document.getElementById("control_AA").value);
    var controls_Aa = parseFloat(document.getElementById("control_Aa").value);
    var controls_aa = parseFloat(document.getElementById("control_aa").value);

    var total_AA = cases_AA + controls_AA;
    var total_Aa = cases_Aa + controls_Aa;
    var total_aa = cases_aa + controls_aa;
    var total_cases = cases_AA + cases_Aa + cases_aa;
    var total_controls = controls_AA + controls_Aa + controls_aa;
    var total = total_cases + total_controls;

    document.getElementById("total_case").value = total_cases;
    document.getElementById("total_control").value = total_controls;
    document.getElementById("total_AA").value = total_AA;
    document.getElementById("total_Aa").value = total_Aa;
    document.getElementById("total_aa").value = total_aa;
    document.getElementById("total").value = total;

    var p_cases = (2 * cases_AA + cases_Aa) / (2 * total_cases);
    var p_controls = (2 * controls_AA + controls_Aa) / (2 * total_controls);
    var p_hat = (2 * total_AA + total_Aa) / (2 * total);
    var Z = 2 * Math.sqrt(total_cases * total_controls) * (p_cases - p_controls) / Math.sqrt(4 * total_AA + total_Aa - 4 * total * p_hat * p_hat);
    var Z_Squared = Math.pow(Z, 2);


    var para = document.getElementById("result");
    para.innerText = "Z squared: " + Math.round(Z_Squared*10000)/10000;
    var result_chi = document.getElementById("result_chi");
    result_chi.innerText = "P-value: " +  Math.round(compute(Z_Squared, 1)*10000)/10000;
    var info = document.getElementById("info");
    info.innerHTML = "P-value calculated based on the &Chi;<sup>2</sup> distribution with 1 degree of freedom";

}