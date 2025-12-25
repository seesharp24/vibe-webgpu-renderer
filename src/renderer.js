import { Vector3 } from './math.js';

export class Renderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.device = null;
        this.context = null;
        this.frameIndex = 0;
        this.raysPerFrame = 5;
        this.totalSamples = 0;

        // Pipelines and resources
        this.computePipeline = null;
        this.blitPipeline = null;
        this.uniformBuffer = null;

        // Ping-pong resources
        this.storageTextureA = null;
        this.storageTextureB = null;
        this.computeBindGroupA = null;
        this.computeBindGroupB = null;
        this.blitBindGroup = null;
    }

    async init() {
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) throw new Error("No GPUAdapter found.");

        this.device = await adapter.requestDevice();
        this.context = this.canvas.getContext('webgpu');
        this.presentationFormat = navigator.gpu.getPreferredCanvasFormat();

        this.context.configure({
            device: this.device,
            format: this.presentationFormat,
            alphaMode: 'premultiplied',
        });

        // Load Shader
        const response = await fetch('src/shaders/tracer.wgsl');
        const shaderCode = await response.text();
        const shaderModule = this.device.createShaderModule({ code: shaderCode });

        // Explicit Layout for Compute Pipeline to handle unfilterable-float
        const computeBindGroupLayout = this.device.createBindGroupLayout({
            entries: [
                {
                    binding: 0,
                    visibility: GPUShaderStage.COMPUTE,
                    buffer: { type: 'uniform' },
                },
                {
                    binding: 1,
                    visibility: GPUShaderStage.COMPUTE,
                    texture: { sampleType: 'unfilterable-float', viewDimension: '2d' }, // History (Read)
                },
                {
                    binding: 2,
                    visibility: GPUShaderStage.COMPUTE,
                    storageTexture: { access: 'write-only', format: 'rgba32float', viewDimension: '2d' }, // Output (Write)
                },
            ],
        });

        const computePipelineLayout = this.device.createPipelineLayout({
            bindGroupLayouts: [computeBindGroupLayout],
        });

        this.computePipeline = this.device.createComputePipeline({
            layout: computePipelineLayout,
            compute: { module: shaderModule, entryPoint: 'main' },
        });

        // Blit Shader: must use textureLoad because unfilterable
        const blitShaderCode = `
            @vertex
            fn vs_main(@builtin(vertex_index) vertexIndex : u32) -> @builtin(position) vec4f {
                var pos = array<vec2f, 3>(
                    vec2f(-1.0, -1.0),
                    vec2f(3.0, -1.0),
                    vec2f(-1.0, 3.0)
                );
                return vec4f(pos[vertexIndex], 0.0, 1.0);
            }

            @group(0) @binding(0) var myTexture : texture_2d<f32>;

            @fragment
            fn fs_main(@builtin(position) FragCoord : vec4f) -> @location(0) vec4f {
                let coords = vec2i(floor(FragCoord.xy));
                let color = textureLoad(myTexture, coords, 0);
                let corrected = pow(color.rgb, vec3f(1.0/2.2));
                return vec4f(corrected, color.a);
            }
        `;
        const blitModule = this.device.createShaderModule({ code: blitShaderCode });

        // Explicit Layout for Blit
        const blitBindGroupLayout = this.device.createBindGroupLayout({
            entries: [
                {
                    binding: 0,
                    visibility: GPUShaderStage.FRAGMENT,
                    texture: { sampleType: 'unfilterable-float', viewDimension: '2d' },
                },
            ]
        });

        const blitPipelineLayout = this.device.createPipelineLayout({
            bindGroupLayouts: [blitBindGroupLayout],
        });

        this.blitPipeline = this.device.createRenderPipeline({
            layout: blitPipelineLayout,
            vertex: { module: blitModule, entryPoint: 'vs_main' },
            fragment: {
                module: blitModule,
                entryPoint: 'fs_main',
                targets: [{ format: this.presentationFormat }]
            },
            primitive: { topology: 'triangle-list' },
        });

        // Uniform Buffer
        this.uniformBufferSize = 256;
        this.uniformBuffer = this.device.createBuffer({
            size: this.uniformBufferSize,
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
        });

        console.log("WebGPU Initialized");

        const observer = new ResizeObserver(entries => {
            for (const entry of entries) {
                const width = entry.contentBoxSize[0].inlineSize;
                const height = entry.contentBoxSize[0].blockSize;
                // Force size to be non-zero
                if (width > 0 && height > 0) {
                    const dpr = window.devicePixelRatio || 1;
                    const displayWidth = Math.floor(width * dpr);
                    const displayHeight = Math.floor(height * dpr);

                    this.canvas.width = Math.max(1, Math.min(displayWidth, this.device.limits.maxTextureDimension2D));
                    this.canvas.height = Math.max(1, Math.min(displayHeight, this.device.limits.maxTextureDimension2D));
                    this.resetAccumulation();
                }
            }
        });
        observer.observe(this.canvas);
    }

    resetAccumulation() {
        if (!this.device || this.canvas.width === 0 || this.canvas.height === 0) return;

        this.frameIndex = 0;
        this.totalSamples = 0;

        if (this.storageTextureA) this.storageTextureA.destroy();
        if (this.storageTextureB) this.storageTextureB.destroy();

        const textureDesc = {
            size: [this.canvas.width, this.canvas.height],
            format: 'rgba32float',
            usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC | GPUTextureUsage.COPY_DST,
        };

        this.storageTextureA = this.device.createTexture(textureDesc);
        this.storageTextureB = this.device.createTexture(textureDesc);

        const computeLayout = this.computePipeline.getBindGroupLayout(0);

        // Bind Group A: Output to A, Read from B
        this.computeBindGroupA = this.device.createBindGroup({
            layout: computeLayout,
            entries: [
                { binding: 0, resource: { buffer: this.uniformBuffer } },
                { binding: 1, resource: this.storageTextureB.createView() }, // History
                { binding: 2, resource: this.storageTextureA.createView() }, // Output
            ],
        });

        // Bind Group B: Output to B, Read from A
        this.computeBindGroupB = this.device.createBindGroup({
            layout: computeLayout,
            entries: [
                { binding: 0, resource: { buffer: this.uniformBuffer } },
                { binding: 1, resource: this.storageTextureA.createView() }, // History
                { binding: 2, resource: this.storageTextureB.createView() }, // Output
            ],
        });
    }

    updateUniforms() {
        const uniforms = new Float32Array(this.uniformBufferSize / 4);

        // Setup Camera looking at origin
        // Box center is approx (0, 1, 0)
        const camPos = new Vector3(0, 1.0, 2.4);
        const target = new Vector3(0, 1.0, 0);
        const up = new Vector3(0, 1, 0);

        let fwd = Vector3.sub(target, camPos);
        fwd = Vector3.normalize(fwd);

        let right = Vector3.normalize(Vector3.cross(fwd, up));
        let newUp = Vector3.normalize(Vector3.cross(right, fwd));

        // Camera Position
        uniforms[0] = camPos.x; uniforms[1] = camPos.y; uniforms[2] = camPos.z;
        // Camera Forward
        uniforms[4] = fwd.x; uniforms[5] = fwd.y; uniforms[6] = fwd.z;
        // Camera Right
        uniforms[8] = right.x; uniforms[9] = right.y; uniforms[10] = right.z;
        // Camera Up
        uniforms[12] = newUp.x; uniforms[13] = newUp.y; uniforms[14] = newUp.z;

        // Frame Count (u32)
        // Layout:
        // Pos: 0-12
        // Fwd: 16-28
        // Right: 32-44
        // Up: 48-60
        // frameCount: 60 (align 4) -> Index 15
        // resolution: 64 (align 8) -> Index 16

        const frameIndView = new Uint32Array(uniforms.buffer);
        frameIndView[15] = this.totalSamples;
        // Resolution
        uniforms[16] = this.canvas.width;
        uniforms[17] = this.canvas.height;
        // Rays Per Frame
        frameIndView[18] = this.raysPerFrame;

        this.device.queue.writeBuffer(this.uniformBuffer, 0, uniforms);
    }

    render() {
        if (!this.device || !this.computePipeline || !this.computeBindGroupA) return;

        this.updateUniforms();

        const commandEncoder = this.device.createCommandEncoder();

        // Ping Pong Logic
        const isEven = (this.frameIndex % 2 === 0);
        const currentBindGroup = isEven ? this.computeBindGroupA : this.computeBindGroupB;
        const currentOutputTexture = isEven ? this.storageTextureA : this.storageTextureB;

        // 1. Compute Pass
        const computePass = commandEncoder.beginComputePass();
        computePass.setPipeline(this.computePipeline);
        computePass.setBindGroup(0, currentBindGroup);
        computePass.dispatchWorkgroups(
            Math.ceil(this.canvas.width / 8),
            Math.ceil(this.canvas.height / 8)
        );
        computePass.end();

        // 2. Render Pass
        const blitBG = this.device.createBindGroup({
            layout: this.blitPipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: currentOutputTexture.createView() },
            ],
        });

        const renderPass = commandEncoder.beginRenderPass({
            colorAttachments: [{
                view: this.context.getCurrentTexture().createView(),
                clearValue: { r: 0, g: 0, b: 0, a: 1 },
                loadOp: 'clear',
                storeOp: 'store',
            }]
        });
        renderPass.setPipeline(this.blitPipeline);
        renderPass.setBindGroup(0, blitBG);
        renderPass.draw(3);
        renderPass.end();

        this.device.queue.submit([commandEncoder.finish()]);

        // Debug Log
        if (this.frameIndex % 60 === 0) {
            console.log("Rendering Frame: ", this.frameIndex, "Total Samples:", this.totalSamples, "Size:", this.canvas.width, this.canvas.height);
        }

        this.frameIndex++;
        this.totalSamples += this.raysPerFrame;
    }
}
