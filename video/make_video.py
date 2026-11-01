"""Render the LinkedIn video (1080x1350, 30 fps, ~43 s) frame by frame and encode it with ffmpeg.
usage: python3 make_video.py viewer_data.json od_graph.json OUT.mp4 [first_frame last_frame]
 viewer_data       : prep_viewer_data.py output of the slices to show
 od_graph           : build_od_graph.py output of the base slices
Needs: playwright (chromium) and ffmpeg. The optional frame range renders only part of the film (for tests)."""
import sys, os, subprocess, tempfile, shutil
from playwright.sync_api import sync_playwright
d1, gr, out = sys.argv[1:4]; rng = [int(v) for v in sys.argv[4:6]] if len(sys.argv) > 5 else None
here = os.path.dirname(os.path.abspath(__file__)); rd = lambda f: open(f, encoding='utf-8').read()
html = rd(f'{here}/stage_template.html').replace('__DATA__', rd(d1)).replace('__GRAPH__', rd(gr))
tmp = tempfile.mkdtemp(prefix='vid_'); page = f'{tmp}/stage.html'; open(page, 'w', encoding='utf-8').write(html)
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1080, 'height': 1350}); errs = []; pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.goto('file://' + page); pg.wait_for_function('typeof window.renderFrame === "function"'); n = pg.evaluate('TOTALF'); a, z = rng if rng else (0, n)
    for f in range(0, z):                      # frames are drawn in order (particles carry state); frames before 'a' are drawn but not saved
        pg.evaluate(f'renderFrame({f})')
        if f >= a: pg.screenshot(path=f'{tmp}/f{f - a:05d}.png', type='png')
        if f % 100 == 0: print('frame', f, '/', z, flush=True)
    b.close()
if errs: print('page errors:', errs); sys.exit(1)
subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '30', '-i', f'{tmp}/f%05d.png', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '17', '-preset', 'medium', '-movflags', '+faststart', out], check=True)
shutil.rmtree(tmp); print('wrote', out)
