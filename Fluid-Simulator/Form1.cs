using System;
using System.Drawing;
using System.Windows.Forms;

namespace Fluid_Simulator
{
    public partial class Form1 : Form
    {
        Fluid fluid;

        int downsample = 4;

        public Form1()
        {
            InitializeComponent();

            fluid = new Fluid(pictureBox1.Width / downsample, 0, 0, 0.1f);

            Paint += Form1_Paint;

            Timer t = new Timer();
            t.Interval = 1000 / 10;
            t.Tick += (s, e) => Invalidate();
            t.Start();
        }

        private void Form1_Paint(object sender, PaintEventArgs e)
        {
            var g = e.Graphics;
            fluid.Step();

            var image = new Bitmap(pictureBox1.Width / downsample, pictureBox1.Width / downsample);

            for (int x = 0; x < pictureBox1.Width / downsample; x++)
            {
                for (int y = 0; y < pictureBox1.Width / downsample; y++)
                {
                    float f = fluid.density[fluid.Get1DIndex(x, y)];
                    int i = (int)f;
                    image.SetPixel(x, y, Color.FromArgb(i, i, i));
                }
            }

            g.DrawImage(image, new Point(0,0));
            image.Dispose();
        }
    }

    class Fluid
    {
        const int iterations = 4;

        int size;
        float timeStep;
        float diffusion;
        float viscosity;

        float[] oldDensity;
        public float[] density;

        float[] velocityX;
        float[] velocityY;

        float[] previouseVelocityX;
        float[] previouseVelocityY;

        public Fluid(int size, int diffusion, int viscosity, float timeStep)
        {
            this.size = size;
            this.diffusion = diffusion;
            this.viscosity = viscosity;
            this.timeStep = timeStep;

            oldDensity = new float[size * size];
            density = new float[size * size];
            velocityX = new float[size * size];
            velocityY = new float[size * size];
            previouseVelocityX = new float[size * size];
            previouseVelocityY = new float[size * size];
        }

        public void AddDye(int x, int y, float ammount)
        {
            int index = Get1DIndex(x, y);
            density[index] += ammount;
        }

        public void AddVelocity(int x, int y, float ammountX, float ammountY)
        {
            int index = Get1DIndex(x, y);
            velocityX[index] += ammountX;
            velocityY[index] += ammountY;
        }

        public int Get1DIndex(int x, int y)
        {
            return y * size + x;
        }

        void Diffuse(int b, ref float[] x, ref float[] oldX, float ammount, float timeStep)
        {
            float a = timeStep * ammount * (size - 2) * (size - 2);
            LinearSolve(b, ref x, ref oldX, a, 1 + 6 * a);
        }


        void Advect(int b, ref float[] d, ref float[] d0, ref float[] velocX, ref float[] velocY, float timeStep)
        {
            float i0, i1, j0, j1;

            float dtx = timeStep * (size - 2);
            float dty = timeStep * (size - 2);

            float s0, s1, t0, t1;
            float tmp1, tmp2, x, y;

            float sizeFloat = size;
            float ifloat, jfloat;
            int i, j;

            for (j = 1, jfloat = 1; j < size - 1; j++, jfloat++)
            {
                for (i = 1, ifloat = 1; i < size - 1; i++, ifloat++)
                {
                    tmp1 = dtx * velocX[Get1DIndex(i, j)];
                    tmp2 = dty * velocY[Get1DIndex(i, j)];
                    x = ifloat - tmp1;
                    y = jfloat - tmp2;

                    if (x < 0.5f) x = 0.5f;
                    if (x > sizeFloat + 0.5f) x = sizeFloat + 0.5f;
                    i0 = (float)Math.Floor(x);
                    i1 = i0 + 1.0f;
                    if (y < 0.5f) y = 0.5f;
                    if (y > sizeFloat + 0.5f) y = sizeFloat + 0.5f;
                    j0 = (float)Math.Floor(y);
                    j1 = j0 + 1.0f;

                    s1 = x - i0;
                    s0 = 1.0f - s1;
                    t1 = y - j0;
                    t0 = 1.0f - t1;

                    int i0i = (int)i0;
                    int i1i = (int)i1;
                    int j0i = (int)j0;
                    int j1i = (int)j1;

                    d[Get1DIndex(i, j)] =
                        s0 * (t0 * d0[Get1DIndex(i0i, j0i)]) + t1 * d0[Get1DIndex(i0i, j1i)] +
                        s1 * (t0 * d0[Get1DIndex(i1i, j0i)]) + t1 * d0[Get1DIndex(i1i, j1i)];
                }
            }

            SetBounds(b, ref d);
        }

        void Project(ref float[] velocX, ref float[] velocY, ref float[] p, ref float[] div)
        {
            for (int j = 1; j < size - 1; j++)
            {
                for (int i = 1; i < size - 1; i++)
                {
                    div[Get1DIndex(i, j)] = -0.5f * (
                                velocX[Get1DIndex(i + 1, j)]
                            - velocX[Get1DIndex(i - 1, j)]
                            + velocY[Get1DIndex(i, j + 1)]
                            - velocY[Get1DIndex(i, j - 1)]
                        ) / size;
                    p[Get1DIndex(i, j)] = 0;
                }
            }

            SetBounds(0, ref div);
            SetBounds(0, ref p);
            LinearSolve(0, ref p, ref div, 1, 6);

            for (int j = 1; j < size - 1; j++)
            {
                for (int i = 1; i < size - 1; i++)
                {
                    velocX[Get1DIndex(i, j)] -=
                        0.5f * (p[Get1DIndex(i + 1, j)]
                        - p[Get1DIndex(i - 1, j)]) * size;
                    velocY[Get1DIndex(i, j)] -=
                        0.5f * (p[Get1DIndex(i, j + 1)]
                        - p[Get1DIndex(i, j - 1)]) * size;
                }
            }

            SetBounds(1, ref velocX);
            SetBounds(2, ref velocY);
        }

        void LinearSolve(int b, ref float[] x, ref float[] oldX, float ammount, float c)
        {
            float cRecip = 1.0f / c;
            for (int k = 0; k < iterations; k++)
            {
                for (int j = 1; j < size - 1; j++)
                {
                    for (int i = 1; i < size - 1; i++)
                    {
                        x[Get1DIndex(i, j)] =
                            (oldX[Get1DIndex(i, j)]
                                + ammount * (x[Get1DIndex(i + 1, j)]
                                        + x[Get1DIndex(i - 1, j)]
                                        + x[Get1DIndex(i, j + 1)]
                                        + x[Get1DIndex(i, j - 1)]
                                        + x[Get1DIndex(i, j)]
                                        + x[Get1DIndex(i, j)]
                                )) * cRecip;
                    }
                }

                SetBounds(b, ref x);
            }
        }

        void SetBounds(int b, ref float[] x)
        {
            for (int k = 1; k < size - 1; k++)
            {
                for (int i = 1; i < size - 1; i++)
                {
                    x[Get1DIndex(i, 0)] = b == 2 ? -x[Get1DIndex(i, 1)] : x[Get1DIndex(i, 1)];
                    x[Get1DIndex(i, size - 1)] = b == 2 ? -x[Get1DIndex(i, size - 2)] : x[Get1DIndex(i, size - 2)];
                }
            }
            for (int k = 1; k < size - 1; k++)
            {
                for (int j = 1; j < size - 1; j++)
                {
                    x[Get1DIndex(0, j)] = b == 1 ? -x[Get1DIndex(1, j)] : x[Get1DIndex(1, j)];
                    x[Get1DIndex(size - 1, j)] = b == 1 ? -x[Get1DIndex(size - 2, j)] : x[Get1DIndex(size - 2, j)];
                }
            }

            x[Get1DIndex(0, 0)] =
                0.33f * (x[Get1DIndex(1, 0)]
                + x[Get1DIndex(0, 1)]
                + x[Get1DIndex(0, 0)]);
            x[Get1DIndex(0, size - 1)] =
                0.33f * (x[Get1DIndex(1, size - 1)]
                + x[Get1DIndex(0, size - 2)]
                + x[Get1DIndex(0, size - 1)]);
            x[Get1DIndex(size - 1, 0)] =
                0.33f * (x[Get1DIndex(size - 2, 0)]
                + x[Get1DIndex(size - 1, 1)]
                + x[Get1DIndex(size - 1, 0)]);
            x[Get1DIndex(size - 1, size - 1)] =
                0.33f * (x[Get1DIndex(size - 2, size - 1)]
                + x[Get1DIndex(size - 1, size - 2)]
                + x[Get1DIndex(size - 1, size - 1)]);
        }

        public void Step()
        {
            Diffuse(1, ref previouseVelocityX, ref velocityX, viscosity, timeStep);
            Diffuse(2, ref previouseVelocityY, ref velocityY, viscosity, timeStep);

            Project(ref previouseVelocityX, ref previouseVelocityY, ref velocityX, ref velocityY);

            Advect(1, ref velocityX, ref previouseVelocityX, ref previouseVelocityX, ref previouseVelocityY, timeStep);
            Advect(2, ref velocityY, ref previouseVelocityY, ref previouseVelocityX, ref previouseVelocityY, timeStep);

            Project(ref velocityX, ref velocityY, ref previouseVelocityX, ref previouseVelocityY);

            Diffuse(0, ref oldDensity, ref density, diffusion, timeStep);
            Advect(0, ref density, ref oldDensity, ref velocityX, ref velocityY, timeStep);
        }
    }
}
