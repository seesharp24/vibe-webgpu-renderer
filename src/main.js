import { Renderer } from './renderer.js';

async function init() {
    if (!navigator.gpu) {
        alert("WebGPU not supported on this browser.");
        return;
    }

    const canvas = document.getElementById('canvas');
    const renderer = new Renderer(canvas);
    const statsDiv = document.getElementById('stats');
    let lastTime = performance.now();
    let frameCount = 0;
    let lastFpsUpdate = 0;

    try {
        await renderer.init();

        function frame() {
            const now = performance.now();
            const dt = now - lastTime;
            lastTime = now;

            renderer.render();

            frameCount++;
            if (now - lastFpsUpdate >= 500) {
                const fps = Math.round(1000 / dt);
                const totalSamples = renderer.totalSamples;
                statsDiv.innerHTML = `
                    FPS: ${fps}<br>
                    Frame: ${renderer.frameIndex}<br>
                    Samples: ${totalSamples}
                `;
                lastFpsUpdate = now;
            }

            requestAnimationFrame(frame);
        }

        requestAnimationFrame(frame);
    } catch (err) {
        console.error("Failed to initialize WebGPU:", err);
    }
}

init();
