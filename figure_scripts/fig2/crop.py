from PIL import Image, ImageDraw, ImageFont
import os

def crop_center(image, crop_width, crop_height):
    width, height = image.size
    left = (width - crop_width) // 2
    top = (height - crop_height) // 2
    right = (width + crop_width) // 2
    bottom = (height + crop_height) // 2
    return image.crop((left, top, right, bottom))

def merge_images(input_folder, output_path, crop_width, crop_height, rows, cols):
    images = []
    titles = ['(a) Torus+Box','(b) Torus-Box','(c) Torus*Box']

    for filename in os.listdir(input_folder):
        if filename.endswith('.png'):
            image_path = os.path.join(input_folder, filename)
            img = Image.open(image_path)
            img = crop_center(img, crop_width, crop_height)
            images.append(img)

    new_image = Image.new('RGB', (cols * crop_width, rows * crop_height+50), color='white')

    font = ImageFont.truetype("arial.ttf", 80)

    draw = ImageDraw.Draw(new_image)
    for i in range(rows):
        for j in range(cols):
            index = i * cols + j
            if index < len(images):
                x_offset = j * crop_width
                y_offset = i * crop_height + 50
                new_image.paste(images[index], (x_offset, y_offset))

                
                title = titles[index]
                _,_,title_width, title_height = draw.textbbox((0, y_offset), title, font=font)
                print(title_width, title_height)
                title_x = x_offset + (crop_width - title_width) // 2
                draw.text((title_x, i * crop_height), title, fill="black", font=font)

    
    new_image.save(output_path)

input_folder = "."

output_path = "fig2.png"

crop_width = 1000
crop_height = 800

rows = 1
cols = 3

merge_images(input_folder, output_path, crop_width, crop_height, rows, cols)
