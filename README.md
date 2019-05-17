# Detection-and-Classification-of-Speed-Limit-Signs

## Introduction.
Traffic sign detection and interpretation is a fundamental aspect of advanced driver assistance
systems, and is essential in the development of fully autonomous vehicles. The objective of this
mini-project was to develop a fully automated system to detect speed limit signs, by applying a
series of imaging techniques to a database of real traffic images. The overall functionality of the
system can be divided into two main components:

1. Detection of a traffic sign within the scene; this represents a “region of interest”
within which other processing needs to be carried out.
2. Classification; Analyse the content of the speed limit sign, and classify the speed
limit.

However sign detection and classification is not merely a trivial task due to adverse features
such as weather, motion, viewpoint variation, physical damage, etc. Each of these adversal
features were encapsulated in the dataset provided and appropriate solutions are detailed
below.

### Preprocessing.
Motion blur causes streaks to appear in the image caused by the movement during the
recording of an image. This effect, has a negative impact on the classification of the images due
to the reduction in clarity of the image, as the high frequency components of the image are lost.
This problem was combated by using a wiener filter on any image that exceeded a threshold
value. This threshold was created using a blur metric based upon the paper “The Blur Effect:”
This algorithm quantifies the amount of blur, by comparing the variation in the pixel intensities of
neighboring pixels in the image and a low pass filtered version of the image.

**Weather** the time of day greatly affected the intensity and contrast of the images. Some images
suffered from a lack of contrast due to night time conditions and others suffered from the glare
of the sun. In order improve the quality of the dark images the mean pixel value was calculated,
if the image in question exceeded this threshold the contrast was adjusted using the imadjust()
function in Matlab. To remove glare or light reflection from the image, a threshold was applied to
the grayscale image only the glare or highest intensity pixels were left in the image. These
pixels intensity in this mask were then scaled back to their neighbors’ values in the actual
image. 

**View Point** each image in the database was taken as a road user approached the sign
therefore the images in the development dataset varied in size. Thus in order to correctly
classify the development dataset, every image was resized to a standardized dimension of
[100x100] pixels.

**Physical Damage** this described any sign that contained imperfections on its face that impeded
the translation of information. Such as dirt or holes that appeared upon the final binary image. In
order to remove these objects the bwareaopen () function was used to morphologically remove
any objects below a certain dimension in the image.

## Sign Detection
Traffic signs are designed to stand out from their natural background, in order for road users to
distinguish their contents. For this reason each sign is positioned at a specific height at the side
of the road, with its own particular shape and color. For example, each speed limit sign has a
circular shape with a strong black number surrounded by a bright red circle.
The detection system created was predicated upon these distinctive features, with the red circle
outlining the signs position within the image. The algorithm begins by converting the image into
the HSV format which separates the image intensity from the Chroma or color information. This
leads to a more accurate threshold of the colour components in the image especially when
dealing with varying intensities across the sign such as shadows.
