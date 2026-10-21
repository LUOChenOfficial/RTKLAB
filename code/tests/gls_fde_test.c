/* Exercise actual static transaction helpers without adding a public API. */
#include "../rtklibdemo5/src/rtkpos.c"

int main(void)
{
    rtk_t *r=(rtk_t *)calloc(1,sizeof(*r)),*full;
    double *liveQ,*liveF;
    rtkim_glscand_t c[3]={{2,{1,2},9.0},{1,{3,0},2.0},{1,{4,0},5.0}};
    if (!r) return 1;
    r->nx=5;r->na=3;
    r->x=zeros(5,1);r->P=zeros(5,5);r->xa=zeros(3,1);r->Pa=zeros(3,3);
    r->x[3]=7.0;r->P[18]=2.0;r->xa[0]=3.0;r->Pa[0]=4.0;
    r->ssat[9].outc[0]=0;r->ssat[9].lock[0]=25;r->ssat[9].vsat[0]=1;
    r->nfix=8;r->com_bias=0.5;r->sol.stat=SOLQ_FIX;r->sol.rr[0]=100.0;
    r->intg.glsQ=zeros(2,2);r->intg.glsF=zeros(3,2);
    r->intg.glsalert=1;r->intg.glsexp=r->intg.glsvalid=12;
    r->intg.pl.hpl=r->intg.pl.vpl=INFINITY;
    r->inthpl=r->intvpl=INFINITY;
    full=glsfdesave(r);
    if (!full) return 2;
    r->x[3]=0;r->P[18]=0;r->xa[0]=0;r->Pa[0]=0;
    r->ssat[9].outc[0]=3;r->ssat[9].lock[0]=0;r->ssat[9].vsat[0]=0;
    r->nfix=0;r->com_bias=2;r->sol.stat=SOLQ_NONE;r->sol.rr[0]=200;
    /* Emulate retry scratch reallocation; snapshot pointers are now stale. */
    free(r->intg.glsQ);free(r->intg.glsF);
    liveQ=r->intg.glsQ=zeros(3,3);liveF=r->intg.glsF=zeros(3,3);
    r->intg.glsexp=50;r->intg.glsvalid=49;
    glsfderollback(r,full);
    if (r->x[3]!=7||r->P[18]!=2||r->xa[0]!=3||r->Pa[0]!=4) return 3;
    if (r->ssat[9].outc[0]!=0||r->ssat[9].lock[0]!=25||!r->ssat[9].vsat[0]) return 4;
    if (r->nfix!=8||r->com_bias!=0.5||r->sol.stat!=SOLQ_FIX||r->sol.rr[0]!=100) return 5;
    if (r->intg.glsQ!=liveQ||r->intg.glsF!=liveF||r->intg.glsready||!r->glsassessed) return 6;
    if (r->intg.glsexp!=12||r->intg.glsvalid!=12||!r->intg.glsalert||!isinf(r->inthpl)) return 7;
    /* No trial may pass with a missing model, any alarm, or nonfinite PL. */
    r->intg.glsready=1;r->intg.glsalert=0;r->inthpl=r->intvpl=1;
    if (!glsfdepassed(r)) return 8;
    r->intg.glsvalid=11;if (glsfdepassed(r)) return 9;
    r->intg.glsvalid=12;r->intg.glsalert=1;if (glsfdepassed(r)) return 10;
    r->intg.glsalert=0;r->intvpl=INFINITY;if (glsfdepassed(r)) return 11;
    qsort(c,3,sizeof(c[0]),glsfdecmp);
    if (c[0].sat[0]!=4||c[1].sat[0]!=3||c[2].order!=2) return 12;
    free(full->x);free(full);
    rtkimfree(r);free(r->x);free(r->P);free(r->xa);free(r->Pa);free(r);
    puts("PASS: failed-history rollback, scratch ownership, unavailable guard, ordered candidate search");
    return 0;
}
