double sum = 0;

#pragma omp parallel for shared(sum, a) reduction(+: sum)
for (int i = 0; i < 10000; i++)
{
    sum += a[i]
}
printf("sum %f",sum);

