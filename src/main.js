import { Renderer } from './renderer.js';

async function init() {
    if (!navigator.gpu) {
        alert("WebGPU not supported on this browser.");
        return;
    }

    const canvas = document.getElementById('canvas');
    const renderer = new Renderer(canvas);

    try {
        await renderer.init();
        
        function frame() {
            renderer.render();
            requestAnimationFrame(frame);
        }
        
        requestAnimationFrame(frame);
    } catch (err) {
        console.error("Failed to initialize WebGPU:", err);
    }
}

init();
