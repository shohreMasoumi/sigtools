BEGIN{OFS="\t"}
{
   split($4,av,",");
   sum=0;
   for(i=0;i<200;i++) sum=sum+av[i];
   print $1,$2,$3,sum/200;
}