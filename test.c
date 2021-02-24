#include<stdio.h>
int main()
{
#if (!defined(__STDC__))
 printf("not_stardard");
#elif defined(__STDC_VERSION__)
 printf("stardard", __STDC_VERSION__);
#else
 printf("old");
#endif
 getchar();
 return 0;
}
