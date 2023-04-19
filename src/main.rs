#[derive(Debug)]
pub struct MatIV {
    data: [f32; 16],
}

impl MatIV {
    pub fn create() -> Self {
        Self {
            data: [
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    pub fn identity(&mut self) {
        self.data = [
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ];
    }

    pub fn multiply(&mut self, mat1: &MatIV, mat2: &MatIV) {
        let mut result = [0.0; 16];

        for i in 0..4 {
            for j in 0..4 {
                let mut element = 0.0;
                for k in 0..4 {
                    element += mat1.data[i * 4 + k] * mat2.data[k * 4 + j];
                }
                result[i * 4 + j] = element;
            }
        }

        self.data = result;
    }

    pub fn scale(&mut self, mat: &MatIV, vec: &[f32; 3]) {
        let [x, y, z] = *vec;

        for i in 0..4 {
            self.data[i * 4] = mat.data[i * 4] * x;
            self.data[i * 4 + 1] = mat.data[i * 4 + 1] * y;
            self.data[i * 4 + 2] = mat.data[i * 4 + 2] * z;
            self.data[i * 4 + 3] = mat.data[i * 4 + 3];
        }
    }

    pub fn translate(&mut self, mat: &MatIV, vec: &[f32; 3]) {
        let [x, y, z] = *vec;

        for i in 0..4 {
            self.data[i * 4] = mat.data[i * 4];
            self.data[i * 4 + 1] = mat.data[i * 4 + 1];
            self.data[i * 4 + 2] = mat.data[i * 4 + 2];
            self.data[i * 4 + 3] = mat.data[i * 4] * x
                + mat.data[i * 4 + 1] * y
                + mat.data[i * 4 + 2] * z
                + mat.data[i * 4 + 3];
        }
    }

    pub fn rotate(&mut self, mat: &MatIV, angle: f32, axis: &[f32; 3]) {
        let [x, y, z] = *axis;
        let s = (angle * std::f32::consts::PI / 180.0).sin();
        let c = (angle * std::f32::consts::PI / 180.0).cos();
        let len = (x * x + y * y + z * z).sqrt();

        let x = x / len;
        let y = y / len;
        let z = z / len;

        let nc = 1.0 - c;

        let a00 = mat.data[0];
        let a01 = mat.data[1];
        let a02 = mat.data[2];
        let a03 = mat.data[3];
        let a10 = mat.data[4];
        let a11 = mat.data[5];
        let a12 = mat.data[6];
        let a13 = mat.data[7];
        let a20 = mat.data[8];
        let a21 = mat.data[9];
        let a22 = mat.data[10];
        let a23 = mat.data[11];

        let b00 = x * x * nc + c;
        let b01 = y * x * nc + z * s;
        let b02 = z * x * nc - y * s;
        let b10 = x * y * nc - z * s;
        let b11 = y * y * nc + c;
        let b12 = z * y * nc + x * s;
        let b20 = x * z * nc + y * s;
        let b21 = y * z * nc - x * s;
        let b22 = z * z * nc + c;

        self.data[0] = a00 * b00 + a10 * b01 + a20 * b02;
        self.data[1] = a01 * b00 + a11 * b01 + a21 * b02;
        self.data[2] = a02 * b00 + a12 * b01 + a22 * b02;
        self.data[3] = a03 * b00 + a13 * b01 + a23 * b02;
        self.data[4] = a00 * b10 + a10 * b11 + a20 * b12;
        self.data[5] = a01 * b10 + a11 * b11 + a21 * b12;
        self.data[6] = a02 * b10 + a12 * b11 + a22 * b12;
        self.data[7] = a03 * b10 + a13 * b11 + a23 * b12;
        self.data[8] = a00 * b20 + a10 * b21 + a20 * b22;
        self.data[9] = a01 * b20 + a11 * b21 + a21 * b22;
        self.data[10] = a02 * b20 + a12 * b21 + a22 * b22;
        self.data[11] = a03 * b20 + a13 * b21 + a23 * b22;
        self.data[12] = mat.data[12];
        self.data[13] = mat.data[13];
        self.data[14] = mat.data[14];
        self.data[15] = mat.data[15];
    }

    pub fn look_at(&mut self, eye: Vec<f32>, center: Vec<f32>, up: Vec<f32>) {
        let [ex, ey, ez] = [eye[0], eye[1], eye[2]];
        let [cx, cy, cz] = [center[0], center[1], center[2]];
        let [ux, uy, uz] = [up[0], up[1], up[2]];

        let fx = cx - ex;
        let fy = cy - ey;
        let fz = cz - ez;

        let rlf = 1.0 / (fx * fx + fy * fy + fz * fz).sqrt();
        let fx = fx * rlf;
        let fy = fy * rlf;
        let fz = fz * rlf;

        let sx = fy * uz - fz * uy;
        let sy = fz * ux - fx * uz;
        let sz = fx * uy - fy * ux;

        let rls = 1.0 / (sx * sx + sy * sy + sz * sz).sqrt();
        let sx = sx * rls;
        let sy = sy * rls;
        let sz = sz * rls;

        let ux = sy * fz - sz * fy;
        let uy = sz * fx - sx * fz;
        let uz = sx * fy - sy * fx;

        self.data[0] = sx;
        self.data[1] = ux;
        self.data[2] = -fx;
        self.data[3] = 0.0;
        self.data[4] = sy;
        self.data[5] = uy;
        self.data[6] = -fy;
        self.data[7] = 0.0;
        self.data[8] = sz;
        self.data[9] = uz;
        self.data[10] = -fz;
        self.data[11] = 0.0;
        self.data[12] = -(sx * ex + sy * ey + sz * ez);
        self.data[13] = -(ux * ex + uy * ey + uz * ez);
        self.data[14] = fx * ex + fy * ey + fz * ez;
        self.data[15] = 1.0;
    }

    pub fn perspective(&mut self, fovy: f32, aspect: f32, near: f32, far: f32) {
        let f = 1.0 / (fovy * std::f32::consts::PI / 360.0).tan();
        let nf = 1.0 / (near - far);

        self.data[0] = f / aspect;
        self.data[1] = 0.0;
        self.data[2] = 0.0;
        self.data[3] = 0.0;
        self.data[4] = 0.0;
        self.data[5] = f;
        self.data[6] = 0.0;
        self.data[7] = 0.0;
        self.data[8] = 0.0;
        self.data[9] = 0.0;
        self.data[10] = (far + near) * nf;
        self.data[11] = -1.0;
        self.data[12] = 0.0;
        self.data[13] = 0.0;
        self.data[14] = 2.0 * far * near * nf;
        self.data[15] = 0.0;
    }

    pub fn transpose(&mut self, mat: &MatIV) {
        let mut result = [0.0; 16];
        for i in 0..4 {
            for j in 0..4 {
                result[i * 4 + j] = mat.data[j * 4 + i];
            }
        }
        self.data = result;
    }

    pub fn inverse(&mut self, mat: &MatIV) -> bool {
        let mut temp = [0.0; 16];

        for i in 0..4 {
            temp[i] = mat.data[i * 4];
            temp[i + 4] = mat.data[i * 4 + 1];
            temp[i + 8] = mat.data[i * 4 + 2];
            temp[i + 12] = mat.data[i * 4 + 3];
        }

        for i in 0..4 {
            let mut t = temp[i * 4];
            if t == 0.0 {
                for j in (i + 1)..4 {
                    t = temp[j * 4 + i];
                    if t != 0.0 {
                        for k in 0..4 {
                            temp[i * 4 + k] += temp[j * 4 + k];
                            temp[j * 4 + k] = temp[i * 4 + k] - temp[j * 4 + k];
                            temp[i * 4 + k] -= temp[j * 4 + k];
                        }
                        break;
                    }
                }
            }
            if t == 0.0 {
                return false;
            }
            for j in 0..4 {
                if j != i {
                    temp[j * 4 + i] /= -t;
                }
            }
            for j in 0..4 {
                for k in (i + 1)..4 {
                    temp[j * 4 + k] += temp[j * 4 + i] * temp[i * 4 + k];
                }
            }
            temp[i * 4 + i] = 1.0 / t;
            for j in 0..4 {
                temp[i * 4 + j] *= t;
            }
        }

        self.data = temp;
        true
    }
}

impl MatIV {
    pub fn new() -> MatIV {
        MatIV { data: [0.0; 16] }
    }
}

impl Default for MatIV {
    fn default() -> Self {
        Self::new()
    }
}

fn main() {
    let mut mat = MatIV::new();
    println!("Initial matrix:\n{:?}", mat);

    mat.identity();
    println!("Identity matrix:\n{:?}", mat);

    let mat1 = MatIV::new();
    let mat2 = MatIV::new();
    mat.multiply(&mat1, &mat2);
    println!("Matrix multiplication:\n{:?}", mat);

    let scale = [1.0, 2.0, 3.0];
    mat.scale(&mat1, &scale);
    println!("Scaled matrix:\n{:?}", mat);
}
