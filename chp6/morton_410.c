#define _CRT_SECURE_NO_WARNINGS 1
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef int newc[3];
typedef int ceng[2];

void read(int a[][8]);                             //读入矩形四叉树
int Tobinary(int k);                               //将十进制的行列号k转换为二进制.
int eq(int m, int n);                              //判断m和n是否相等，相等则返回1，否则返回0.
void judge_four(int a[][8], newc b[][8], ceng c[4]);//判断数组c表示的矩形范围内的值是否一样，若一样就更新数组b.
void judge_ten(int a[][8], int b[][8], ceng c[4]);
void checkcombine_four(int a[][8], newc b[][8]);   //将属性值一样的单元进行合并.
void checkcombine_ten(int a[][8], int b[][8]);
void Change(int m[][8], newc n[][8]);              //将四进制的M码按照规定格式输出.
void four_decimal(int b[][8]);                     //四进制编码.
void ten_decimal(int c[][8]);                      //十进制编码.
void output_i(int a[][8]);                         //输出int型数组.
void output_c(newc a[][8]);                        //输出newc型数组.

void main()
{
    int a[8][8], b[8][8], c[8][8];
    newc bb[8][8];
    read(a);
    printf("矩阵四叉树为：\n");
    output_i(a);

    four_decimal(b);
    Change(b, bb);
    checkcombine_four(a, bb);
    printf("\n四叉树对应的四进制编码为:\n");
    output_c(bb);

    ten_decimal(c);
    checkcombine_ten(a, c);
    printf("\n四叉树对应的十进制编码为:\n");
    output_i(c);
}

void read(int a[][8])
{
    FILE* fp = fopen("squad_tree2.txt", "r");
    if (!fp)
        exit;
    else
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                fscanf(fp, "%d", &a[i][j]);
    fclose(fp);
}

int Tobinary(int k)
{
    int s[10], rem, i = 0, t = 0;
    do
    {
        rem = k % 2;
        k = k / 2;
        s[i++] = rem;
    } while (k != 0);   //当十进制数是0时也要进行一遍此循环，所以必须用do……while循环，而不是while循环
    for (int j = --i; j >= 0; j--)
    {
        t += s[j] * pow(10.0, j);
    }
    return t;
}

int eq(int m, int n)
{
    if (m == n)
        return 1;
    else
        return 0;
}

void judge_four(int a[][8], newc b[][8], ceng c[4])
{
    for (int i = 0; i < 4; i++)
    {
        int w = 0;
        for (int m = c[i][0]; m < c[i][0] + 2; m++)
            for (int n = c[i][1]; n < c[i][1] + 1; n++)
                w += eq(a[m][n], a[m][n + 1]);
        if (w == 2)//4个值属性一样
        {
            for (int m = c[i][0]; m < c[i][0] + 2; m++)
                for (int n = c[i][1]; n < c[i][1] + 2; n++)
                {
                    *b[m][n] = *b[(c[i][0])][(c[i][1])];
                    *(b[m][n] + 1) = *(b[(c[i][0])][(c[i][1])] + 1);
                    *(b[m][n] + 2) = 99999999;
                }
        }
    }
}

void judge_ten(int a[][8], int b[][8], ceng c[4])
{
    for (int i = 0; i < 4; i++)
    {
        int w = 0;
        for (int m = c[i][0]; m < c[i][0] + 2; m++)
            for (int n = c[i][1]; n < c[i][1] + 1; n++)
                w += eq(a[m][n], a[m][n + 1]);
        if (w == 2)//4个值属性一样
        {
            for (int m = c[i][0]; m < c[i][0] + 2; m++)
                for (int n = c[i][1]; n < c[i][1] + 2; n++)
                {
                    b[m][n] = b[(c[i][0])][(c[i][1])];
                }
        }
    }
}

void checkcombine_ten(int a[][8], int b[][8])
{
    //第一层
    ceng c[4] = { { 0, 0 }, { 0, 4 }, { 4, 0 }, { 4, 4 } };
    for (int i = 0; i < 4; i++)
    {
        int w = 0;
        for (int m = c[i][0]; m < c[i][0] + 4; m++)
            for (int n = c[i][1]; n < c[i][1] + 3; n++)
                w += eq(a[m][n], a[m][n + 1]);
        if (w == 12)//16个值属性一样
        {
            for (int m = c[i][0]; m < c[i][0] + 4; m++)
                for (int n = c[i][1]; n < c[i][1] + 4; n++)
                {
                    b[m][n] = b[(c[i][0])][(c[i][1])];
                }
        }
    }
    //第二层
    ceng d[4] = { { 0, 0 }, { 0, 2 }, { 2, 0 }, { 2, 2 } },
        e[4] = { { 0, 4 }, { 0, 6 }, { 2, 4 }, { 2, 6 } },
        f[4] = { { 4, 0 }, { 4, 2 }, { 6, 0 }, { 6, 2 } },
        g[4] = { { 4, 4 }, { 4, 6 }, { 6, 4 }, { 6, 6 } };
    judge_ten(a, b, d);
    judge_ten(a, b, e);
    judge_ten(a, b, f);
    judge_ten(a, b, g);
}

void checkcombine_four(int a[][8], newc b[][8])
{
    //第一层
    ceng c[4] = { { 0, 0 }, { 0, 4 }, { 4, 0 }, { 4, 4 } };
    for (int i = 0; i < 4; i++)
    {
        int w = 0;
        for (int m = c[i][0]; m < c[i][0] + 4; m++)
            for (int n = c[i][1]; n < c[i][1] + 3; n++)
                w += eq(a[m][n], a[m][n + 1]);
        if (w == 12)//16个值属性一样
        {
            for (int m = c[i][0]; m < c[i][0] + 4; m++)
                for (int n = c[i][1]; n < c[i][1] + 4; n++)
                {
                    *b[m][n] = *b[(c[i][0])][(c[i][1])];
                    *(b[m][n] + 1) = 99999999;
                    *(b[m][n] + 2) = 99999999;
                }
        }
    }
    //第二层
    ceng d[4] = { { 0, 0 }, { 0, 2 }, { 2, 0 }, { 2, 2 } },
        e[4] = { { 0, 4 }, { 0, 6 }, { 2, 4 }, { 2, 6 } },
        f[4] = { { 4, 0 }, { 4, 2 }, { 6, 0 }, { 6, 2 } },
        g[4] = { { 4, 4 }, { 4, 6 }, { 6, 4 }, { 6, 6 } };
    judge_four(a, b, d);
    judge_four(a, b, e);
    judge_four(a, b, f);
    judge_four(a, b, g);
}

void Change(int m[][8], newc n[][8])
{
    int t[3];
    int q;
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
        {
            q = m[i][j];
            t[0] = q / 100;
            q = q % 100;
            t[1] = q / 10;
            q = q % 10;
            t[2] = q;
            *n[i][j] = *t;             //数组赋值，数组名称不能直接做左值
            *(n[i][j] + 1) = *(t + 1);
            *(n[i][j] + 2) = *(t + 2);
        }
}

void four_decimal(int b[][8])
{
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
        {
            int m = Tobinary(i);
            int n = Tobinary(j);
            b[i][j] = 2 * m + n;
        }
}

void ten_decimal(int c[][8])
{
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
        {
            int m = Tobinary(i);
            int n = Tobinary(j);
            int t[8];
            t[0] = m / 1000;
            m = m % 1000;
            t[2] = m / 100;
            m = m % 100;
            t[4] = m / 10;
            m = m % 10;
            t[6] = m;

            t[1] = n / 1000;
            n = n % 1000;
            t[3] = n / 100;
            n = n % 100;
            t[5] = n / 10;
            n = n % 10;
            t[7] = n;

            int y = 0;
            for (int w = 0; w < 8; w++)
            {
                y += t[w] * pow(2.0, 7 - w);
            }
            c[i][j] = y;
        }
}

void output_i(int a[][8])
{
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            printf("%6d", a[i][j]);
        }
        printf("\n");
    }
}

void output_c(newc a[][8])
{
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
            if (*(a[i][j] + 1) == 99999999 && *(a[i][j] + 2) == 99999999)
                printf("%6d", *a[i][j]);
            else if (*(a[i][j] + 2) == 99999999)
                printf("%5d%d", *a[i][j], *(a[i][j] + 1));
            else
                printf("%4d%d%d", *a[i][j], *(a[i][j] + 1), *(a[i][j] + 2));
        printf("\n");
    }
}
