from django.contrib import admin
from django.urls import include, re_path,path
from django.conf import settings
from django.conf.urls.static import static
from .views import QueryElement,QueryElement_tutorial



urlpatterns = [
    re_path(r'^tutorial1$', QueryElement_tutorial.as_view(), name='tutorial_1'),
    re_path('', QueryElement.as_view(), name='query'),
    ]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)