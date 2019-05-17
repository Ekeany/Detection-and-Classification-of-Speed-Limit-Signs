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

#### Preprocessing.
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

<p align="center">
  <img width="644" height="344" src="/Images_/image1.PNG">
</p>

By limiting the range in HSV space the image was converted to a binary image highlighting the
red Chroma. To isolate the sign in the binary image the imfindcircles() function was then applied
to the image. This function is based upon the circular Hough transform and returns the position
and radius of any circle found in the image. The circle with the largest radius found in the upper
half of the image was taken as the correct location of the sign. This region of interest was then
extracted using basic geometry and the imcrop() function supplied by Matlab. 

<p align="center">
  <img width="780" height="630" src="/Images_/image2.PNG">
</p>

This algorithm alone was only able to extract 50% of the images in the stress dataset as some
signs are surrounded by a square orange background that is also converted to binary resulting
in a square sign. Thus if the imfindcircles() function fails to find any circles in the binary image
the algorithm then searches for squares. Using the regionprops() function the area,
boundingbox , major and minor axis were extracted for every blob in the image.

A square is defined as a shape that has four equal sides at right angles to each other. To
classify a blob as a square the ratio between the minor and major axis of the blob was
calculated along with the ratio of the blob’s pixel area to the area produced by multiplying both
the minor and major axis. If both of these ratios were contained within a range close to one then
the blob was classified as a square object. These square blobs are then filtered based upon
their position relative to the center of the image and the most appropriate square is then.
extracted from the image.

<p align="center">
  <img width="632" height="431" src="/Images_/image3.PNG">
</p>

<p align="center">
  <img width="757" height="322" src="/Images_/image4.PNG">
</p>

### Classification.
Within this mini project classification was conducted by comparing the extracted features of the
processed image in question to five gold standard images. Classification begins by first
selecting a region of interest from each image in the development dataset. Each image was
firstly rescaled to [100 100] pixels and each of the preprocessing steps described above were
implemented when appropriate, to increase contrast.

<p align="center">
  <img width="509" height="290" src="/Images_/image5.PNG">
</p>

Using the imfindcircles() function on the newly preprocessed grayscale images the position and
size of each circle in the image was found. As the sign was the subject of focus from the
development set images. The circle that contained the boundary of the speed sign was
assumed to be at the center of the image with a radius of 15-45 pixels. Therefore the distance
from each circle to the center of the image was calculated and the circle closest to the center
was extracted from the image. 

<p align="center">
  <img width="505" height="270" src="/Images_/image6.PNG">
</p>

A mask was firstly placed around the circle to crop the sign from its natural background.
Then a region of interest was selected from the image using some basic geometry paired with
the imcrop() function to select the left hand semi-circle of the sign. This portion of the image was
selected as it contains the unique digit in each speed sign. This process was also conducted for
each of the images in the gold standard dataset and each image was then rescaled to the exact
same size an converted to binary. 

<p align="center">
  <img width="623" height="257" src="/Images_/image7.PNG">
</p>

<p align="center">
  <img width="629" height="167" src="/Images_/image8.PNG">
</p>

From previous experience the most accurate feature descriptor of an image for the purpose of
object detection is the histogram of oriented gradients. This theory essentially states that a local
object’s appearance and shape within an image can be described by the distribution of intensity
gradients or edge directions. The algorithm begins by dividing the image into small connected
regions called cells. For the pixels within each cell, a histogram of gradient directions is
compiled. This feature vector was created for all the processed images in the development
dataset and compared using the Euclidean distance to each feature vector in the Gold standard
dataset the image with the shortest distance was then chosen to be the class of the unknown
image. 

<p align="center">
  <img width="517" height="423" src="/Images_/image9.PNG">
</p>

## Results.
### Detection.
The algorithm used to detect the signs from the stress dataset achieved an accuracy of 93.75%
on the images in the stress dataset, only failing to detect one sign from 16 images. This failure
is based on the nonexistent objects in the binary image. For example the sign in the third image
was positioned too far in the background of the image to be detected by the HSV range
threshold.

<p align="center">
  <img width="654" height="285" src="/Images_/image10.PNG">
</p>

### Classification.
The original classification accuracy using the minimum Euclidean distance between feature
vectors of an unknown image and an image from the gold standard are detailed below. The
average accuracy over the entire development set was recorded at 53% with the 20kmh signs
achieving the highest accuracy at 70% and the 100kmh signs achieving the lowest accuracy at
39%. This discrepancy in accuracy can be attributed to the difference in the quality of the
images in the dataset and the distinct shape of the digits themselves. The 20kmh signs had the
highest contrast and readability whereas the 100kmh and 80kmh signs suffered greatly from
motion blur and lighting issues. Using the exact same classification algorithm the average
classification accuracy on the stress dataset was only 25%.

<p align="center">
  <img width="649" height="531" src="/Images_/image11.PNG">
</p>

In order to obtain a classification for the 80kmh and 100kmh development folders two new
images had to be downloaded from the internet and resized to their appropriate proportions. It
was noted that by using various different signs caused the accuracy to fluctuate wildly. The
highest accuracy was achieved with an eighty kilometer sign that had tall and narrow digits
which allowed it to be distinguished from the thirty and fifty kilometers signs. To better
understand these fluctuations the hog feature vectors were visualized in vector space using
dimensionality reduction techniques. In doing so it was evident that the gold standard images
were not a fit representation of the development dataset. 

<p align="center">
  <img width="648" height="308" src="/Images_/image12.PNG">
</p>
<p align="center">
  <img width="633" height="504" src="/Images_/image13.PNG">
</p>

Thus by simply taking the average feature vector of selection of training examples from each
folder 5 new points were created that better encapsulated each class. This new method was
tested on 20% of the images in the dataset again using the minimum Euclidean distance from
the unknown feature vector to classify the image. This improved the accuracy of the
development dataset to 66% and drastically improved the classification accuracy of the stress
dataset to 62.5%. The new results are displayed below.

<p align="center">
  <img width="612" height="494" src="/Images_/image14.PNG">
</p>

## Conclusion.
The automatic detection and classification of traffic signs is not a trivial task. In real world
applications many uncontrollable elements can alter the performance of the program. This was
clearly evident in the development dataset where some folders achieved a lower accuracy score
due to their contrast being lowered by motion blur, lighting and sign degradation. A number of
preprocessing steps were implemented however a number of these images were unsalvageable
using the techniques described in this report. Despite this the program achieved respectable
sign detection rates at 92% in the stress dataset however the classification accuracy of the
signs could be improved. The histogram of oriented gradients produces a reliable feature
description of each image. Despite this, merely using the Euclidean distance between golden
standard points in feature space is not the most accurate method for classification. This
accuracy could have been drastically increased by introducing machine learning into program. 

## Reference.
1. https://www.mathworks.com/matlabcentral/fileexchange/28689-hog-descriptor-for-matlab
2. https://www.mathworks.com/matlabcentral/fileexchange/24676-image-blur-metric
3. https://nuigalway.blackboard.com/webapps/blackboard/execute/content/file?cmd=view&content_id=_1687003_1&course_id=_98245_1

