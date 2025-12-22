export class Vector3 {
    constructor(x = 0, y = 0, z = 0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    static normalize(v) {
        const len = Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        if (len > 0) {
            return new Vector3(v.x / len, v.y / len, v.z / len);
        }
        return new Vector3(0, 0, 0);
    }

    static sub(a, b) {
        return new Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    static cross(a, b) {
        return new Vector3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        );
    }

    // Quick helper for buffer data
    toArray() {
        return [this.x, this.y, this.z];
    }
}

// Minimal matrix for Camera
export class Matrix4 {
    constructor() {
        this.elements = new Float32Array(16);
        this.identity();
    }

    identity() {
        this.elements.fill(0);
        this.elements[0] = 1;
        this.elements[5] = 1;
        this.elements[10] = 1;
        this.elements[15] = 1;
        return this;
    }

    // A simple lookAt implementation
    lookAt(eye, center, up) {
        const z = Vector3.normalize(Vector3.sub(eye, center));
        const x = Vector3.normalize(Vector3.cross(up, z));
        const y = Vector3.normalize(Vector3.cross(z, x));

        const te = this.elements;
        te[0] = x.x; te[4] = y.x; te[8] = z.x; te[12] = eye.x;
        te[1] = x.y; te[5] = y.y; te[9] = z.y; te[13] = eye.y;
        te[2] = x.z; te[6] = y.z; te[10] = z.z; te[14] = eye.z;
        te[3] = 0; te[7] = 0; te[11] = 0; te[15] = 1;

        return this;
    }
}
