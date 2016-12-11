void nuc_loop( char* p, int stride, char* q, int u, int v )
{
    u *= stride ;
    v *= stride ;

    while( u < v ) {
        char a = q[ u ] ;
        char b = q[ u + stride ] ;
        char a1 = a ? 0x10 << (a&3) : 0xf0 ;
        char b1 = b ? 0x1  << (b&3) : 0xf  ;
        *p++ = a1 | b1 ;
        u += stride+stride ;
    }
    if( u == v ) {
        char a = q[ u ] ;
        char a1 = a ? 0x10 << (a&3) : 0xf0 ;
        *p = a1 ;
    }
}

void nuc_loop_asc( char* p, int stride, char* q, int u, int v )
{
    u *= stride ;
    v *= stride ;

    while( u <= v ) {
        char a = q[ u ] ;
        *p++ = a == 0 ? 'N' : (a&3) == 0 ? 'A' : (a&3) == 1 ? 'C' : (a&3) == 2 ? 'G' : 'T' ;
        u += stride ;
    }
    *p = 0 ;
}

void qual_loop( char* p, int stride, char* q, int u, int v )
{
    u *= stride ;
    v *= stride ;
    while( u <= v ) {
        *p++ = (q[u] >> 2) & 0x3f ;
        u += stride ;
    }
}

void qual_loop_asc( char* p, int stride, char* q, int u, int v )
{
    u *= stride ;
    v *= stride ;
    while( u <= v ) {
        *p++ = 33 + ((q[u] >> 2) & 0x3f) ;
        u += stride ;
    }
    *p = 0 ;
}

int int_loop( char* p, int x )
{
    *p++ = ':' ;
    if( x == 0 ) {
        *p = '0' ;
        return 2 ;
    }
    char *q = p ;
    while( x > 0 ) {
        *q++ = '0' + x % 10 ;
        x /= 10 ;
    }
    int r = q-p ;
    --q ;
    while( p < q ) {
        char c = *p ;
        *p++ = *q ;
        *q-- = c ;
    }
    return r+1 ;
}


