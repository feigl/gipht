double arc(double a, double b)
{
    double c;
    
    /*
     * given phase gradient values a and b each in radians on [-pi,pi]
     * return arc(a,b) on [0 pi]
     * see eq 2.3.13 page 19, Mardia and Jupp [2000]
     *
     * Kurt Feigl 2010-DEC-03
     */
    
    c = M_PI - abs(M_PI-abs(a-b));
    return(c);
}
