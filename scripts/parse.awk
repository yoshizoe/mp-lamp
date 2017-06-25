/#Database:/ {N=$3; n=$5; N1=$7;}
/1st phase end:/{enumtime=$10;}
/2nd phase end:/{ntest=$6; sigtime=$12;}
/number of significant patterns=/{nsig=$6;}
/3rd phase end:/{sigtesttime=$8}
END{
    printf("%d %d %d %d %d %.2f %.2f\n", N, N1, n, ntest, nsig, enumtime, sigtime + sigtesttime);
}
